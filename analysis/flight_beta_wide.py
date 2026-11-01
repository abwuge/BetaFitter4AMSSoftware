#!/usr/bin/env python3
"""Flight-only mass fits over 2.5--200 GV; reserved RICH evaluation.

Three independently trained rigidity blocks retain the original four-bin
mass-likelihood implementation. No RICH measurement enters extraction or fit.
"""
import argparse
import json
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np
import uproot

import train_flight_mass_beta as tr
import rich_transport as rt
import evaluate_rich_transport as ev

DEST = rt.ROOT / 'results/flight_beta_wide'
CACHE = Path('/tmp/betaEstimator-flight-beta-wide')
BLOCKS = [np.array(x, float) for x in
          [[2.5, 3, 4, 5, 6], [6, 8, 10, 15, 20], [20, 30, 50, 100, 200]]]
R_EDGES = np.unique(np.concatenate(BLOCKS))
BETA_EDGES = np.array([.78, .85, .90, .95, .98, .99, .995, .998, 1.02])


def set_edges(edges):
    rt.rc.R_EDGES = np.array(edges, float)
    rt.rc.N_BINS = len(edges) - 1


def selection(a):
    d = rt.center_features(a)
    q = np.where(a['tk_qin'][:, 0, 2] < 2.5, a['tk_q'][:, 1], a['tk_qin'][:, 0, 2])
    ok = np.isfinite(q) & (q >= 5.5) & (q < 6.5) & (a['tof_bpatt'] == 4444)
    ok &= np.isfinite(d['u']) & (d['u'] > 0) & np.all(np.isfinite(d['continuous']), axis=1)
    ok &= (d['bin'] >= 0) & (d['bin'] < len(R_EDGES) - 1)
    ok &= np.all((a['tof_barid'] >= 0) & (a['tof_barid'] < 10), axis=1)
    ok &= np.all(np.isfinite(a['tk_pos'][:, 1:8]), axis=(1, 2))
    ok &= (a['tk_pos'][:, 4, 2] > 0) & (a['tk_pos'][:, 5, 2] < 0)
    rel = np.abs(a['tk_erigidity1'][:, 2, 1] * d['R'])
    chi = a['tk_chis1'][:, 2, 2, 1, :]
    ok &= np.isfinite(rel) & (rel > 0) & (rel < .20)
    ok &= np.all(np.isfinite(chi) & (chi > 0) & (chi <= 10), axis=1)
    return ok


def extract_part(task):
    source, role, path = task
    set_edges(R_EDGES)
    branches = tr.TRAIN_BRANCHES if role == 'train' else list(dict.fromkeys(ev.MEASURED + ev.SELECTION + ['rich_q']))
    part = CACHE / role / ('%04d.npz' % source)
    note = part.with_suffix('.json')
    signature = dict(path=path, source=source, role=role, branches=branches,
                     script_sha256=ev.digest(Path(__file__)), edges=R_EDGES.tolist())
    if part.exists() and note.exists():
        audit = json.loads(note.read_text())
        if audit['signature'] == signature and audit['complete']:
            return str(part), audit
    pieces, raw = [], 0
    with uproot.open(path) as f:
        for a in f['amstreea'].iterate(branches, library='np', step_size='100 MB'):
            raw += len(a['run'])
            rows = np.flatnonzero(selection(a))
            if role == 'train':
                d = rt.center_features(ev.sub(a, rows))
                d.update(run=a['run'][rows], event=a['event'][rows],
                         source=np.full(len(rows), source, np.int32),
                         relR=np.abs(a['tk_erigidity1'][rows, 2, 1] * d['R']))
                if set(d) != tr.TRAIN_KEYS:
                    raise ValueError('Training fields must contain measurements only')
            else:
                d = ev.sub(a, rows)
                d['source'] = np.full(len(rows), source, np.int32)
            pieces.append(d)
    data = ev.concatenate(pieces)
    temporary = part.with_suffix('.partial.npz')
    np.savez_compressed(str(temporary), **data)
    temporary.replace(part)
    audit = dict(signature=signature, raw=raw, selected=len(data['run']), complete=True)
    ev.dump(note, audit)
    return str(part), audit


def extract(role, workers):
    if role == 'validation':
        check_freeze()
    DEST.mkdir(parents=True, exist_ok=True)
    (CACHE / role).mkdir(parents=True, exist_ok=True)
    tasks = [s for s in tr.sources() if s[1] == role]
    pieces, audit = [], []
    with mp.Pool(workers) as pool:
        for j, (path, record) in enumerate(pool.imap(extract_part, tasks), 1):
            with np.load(path, allow_pickle=False) as f:
                pieces.append({k: f[k] for k in f.files})
            audit.append(record)
            if j % 20 == 0 or j == len(tasks):
                print(role, j, '/', len(tasks), 'selected', sum(a['selected'] for a in audit), flush=True)
    data = ev.concatenate(pieces)
    ids = np.rec.fromarrays([data['run'], data['event']], names='run,event')
    if len(np.unique(ids)) != len(ids):
        raise ValueError('Duplicate event identities')
    other_runs = {int(Path(s[2]).name.split('_')[0]) for s in tr.sources() if s[1] != role}
    if set(data['run'].tolist()) & other_runs:
        raise ValueError('Training/validation runs overlap')
    set_edges(R_EDGES)
    bins = data['bin'] if role == 'train' else rt.center_features(data)['bin']
    np.savez_compressed(str(CACHE / (role + '.npz')), **data)
    ev.dump(DEST / (role + '_extraction.json'), dict(complete=True, n=len(ids),
            raw=sum(a['raw'] for a in audit), sources=audit, edges=R_EDGES,
            bins=np.bincount(bins, minlength=len(R_EDGES)-1), unique_events=True, disjoint_runs=True))


def dependencies():
    return {str(Path(p).resolve().relative_to(rt.ROOT)): ev.digest(p) for p in
            [__file__, tr.__file__, rt.__file__, ev.__file__, rt.rc.__file__, rt.er.__file__,
             rt.ROOT / 'config/data_sources.conf']}


def train():
    if (DEST / 'freeze.json').exists():
        raise RuntimeError('Already frozen; no refit after RICH evaluation')
    with np.load(str(CACHE / 'train.npz'), allow_pickle=False) as f:
        if set(f.files) != tr.TRAIN_KEYS:
            raise ValueError('Unexpected training fields')
        data = {k: f[k] for k in f.files}
    summary = []
    for i, edges in enumerate(BLOCKS):
        set_edges(edges)
        d = ev.sub(data, (data['R'] >= edges[0]) & (data['R'] < edges[-1]))
        d['bin'] = np.searchsorted(edges, d['R'], side='right') - 1
        counts = np.bincount(d['bin'], minlength=4)
        if min(counts) < 100:
            raise RuntimeError('Insufficient training coverage: ' + str(counts))
        print('BLOCK', i, 'edges', edges, 'N', len(d['u']), 'counts', counts, flush=True)
        raw = tr.fit_raw(d)
        print('raw fractions', raw['fraction'], flush=True)
        state, matrix = tr.features(d)
        # Check actual wider-domain likelihood gradient before optimization.
        density = rt.rc.MixtureDensity(d, matrix, raw['fraction'][d['bin']])
        p0 = np.r_[np.zeros(matrix.shape[1]), raw['logsigma']]
        direction = np.random.RandomState(13).normal(size=len(p0))
        direction[0] = 0
        direction /= np.linalg.norm(direction)
        eps = 1e-6
        numeric = (density.value_gradient(p0+eps*direction)[0]-density.value_gradient(p0-eps*direction)[0])/(2*eps)
        analytic = density.value_gradient(p0)[1].dot(direction)
        if abs(numeric-analytic) > 1e-5 * max(1., abs(analytic)):
            raise RuntimeError('Directional gradient mismatch')
        params, record = tr.fit_correction(d, matrix, raw, raw['fraction'])
        record.update(edges=edges, n=len(d['u']), counts=counts, raw_spectrum=raw,
                      gradient_difference=abs(numeric-analytic))
        model = dict(state, parameters=params, rigidity_edges_GV=edges)
        # Replay from serialized state, not just the in-memory training matrix.
        path = DEST / ('block_%d_model.json' % i)
        ev.dump(path, model)
        reload = rt.arrays(json.loads(path.read_text()))
        replay = tr.feature_matrix(d, reload)
        record['serialization_max_difference'] = float(np.max(np.abs(matrix.dot(params[:matrix.shape[1]]) - replay.dot(params[:replay.shape[1]]))))
        if record['serialization_max_difference'] > 1e-12:
            raise RuntimeError('Serialization mismatch')
        summary.append(record)
        ev.dump(DEST / 'training.json', dict(blocks=summary, complete=len(summary)==len(BLOCKS),
                RICH_used=False, MC_used=False, coefficients='independent per rigidity block',
                fraction_status='training-only, response-dependent nuisance; not abundance measurement'))
        print('FIT COMPLETE', i, record['iterations'], record['rms_correction'], flush=True)
    ev.dump(DEST / 'freeze.json', dict(utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
            dependencies=dependencies(), model_sha256={str(i): ev.digest(DEST/('block_%d_model.json'%i)) for i in range(3)},
            train_sha256=ev.digest(CACHE/'train.npz'), material_sha256=ev.digest(rt.MAP),
            rigidity_edges_GV=R_EDGES, beta_edges=BETA_EDGES,
            propagation='C12 expected inverse-beta delay, additive to measured inverse beta; retain measured beta>1',
            reference='same RICH quality plus 4<Q<8; no TOF-RICH agreement or mass cut',
            models='three independent blocks, nominal training-only fractions; no RICH selection'))
    print('FROZEN', flush=True)


def check_freeze():
    freeze = json.loads((DEST/'freeze.json').read_text())
    if freeze['dependencies'] != dependencies():
        raise RuntimeError('Frozen dependency changed')
    for i, digest in freeze['model_sha256'].items():
        if ev.digest(DEST/('block_%s_model.json'%i)) != digest:
            raise RuntimeError('Model changed')
    if freeze['material_sha256'] != ev.digest(rt.MAP):
        raise RuntimeError('Material map changed')
    return freeze


def predict(d):
    corrected = np.full(len(d['u']), np.nan)
    for i, edges in enumerate(BLOCKS):
        set_edges(edges)
        rows = (d['R'] >= edges[0]) & (d['R'] < edges[-1])
        local = ev.sub(d, rows)
        local['bin'] = np.searchsorted(edges, local['R'], side='right') - 1
        state = rt.arrays(json.loads((DEST/('block_%d_model.json'%i)).read_text()))
        matrix = tr.feature_matrix(local, state)
        corrected[rows] = local['u'] - matrix.dot(state['parameters'][:matrix.shape[1]])
    return corrected


def evaluate():
    check_freeze()
    with np.load(str(CACHE/'validation.npz'), allow_pickle=False) as f:
        a = {k:f[k] for k in f.files}
    set_edges(R_EDGES)
    d = rt.center_features(a)
    corrected = predict(d)
    if not np.isfinite(corrected).all():
        raise RuntimeError('Missing covered-domain prediction')
    geometry = rt.MaterialMap().build_columns(a['tk_pos'], a['tk_dir'])
    # Compute material delay from measured rigidity and known mass, rather
    # than assigning an imaginary mass to noisy measurements with beta>1.
    # The same additive u delay is used for old and new estimates.
    delays, valid = [], geometry['valid'].copy()
    for mass in (11.17793, 12.11221):
        nominal_u = np.sqrt(1 + (mass/(6*d['R']))**2)
        physical = rt.propagate(nominal_u, d['R'], geometry, mass_override=mass)
        delay = 1/physical['beta_planes'][:,-1]-nominal_u
        valid &= physical['valid']
        delays.append(delay)
    raw = 1/(d['u']+delays[0])
    new = 1/(corrected+delays[0])
    reference = a['rich_beta'][:,0]
    ring = rt.rich_quality(a)
    charge = (a['rich_npe'][:,2] > 0)  # Q is read separately below if not in old MEASURED allowlist.
    if 'rich_q' not in a:
        raise RuntimeError('Missing independent RICH charge')
    charge &= (a['rich_q'][:,1] > 4) & (a['rich_q'][:,1] < 8)
    mask = ring & charge & valid & np.isfinite(raw) & np.isfinite(new) & (raw>0) & (new>0)
    np.savez_compressed(str(DEST/'predictions.npz'), run=a['run'], event=a['event'], source=a['source'],
                        R=d['R'], rich_beta=reference, raw_u=d['u'], corrected_u=corrected,
                        raw_beta_rich=raw, new_beta_rich=new, reference_mask=mask,
                        delay_u_C12=delays[0], delay_u_C13=delays[1], NaF=a['rich_NaF'])
    rows = {}
    for axis, values, edges in [('rigidity',d['R'],R_EDGES), ('beta',reference,BETA_EDGES)]:
        records = []
        for lo, hi in zip(edges[:-1],edges[1:]):
            select = mask & (values>=lo) & (values<hi)
            n = int(select.sum())
            row = dict(low=lo, high=hi, n=n)
            if n>=2:
                row.update(ev.compare(raw[select],new[select],reference[select]))
                row['x'] = np.median(values[select])
                row['bootstrap'] = ev.bootstrap(raw[select],new[select],reference[select],a['source'][select])
            records.append(row)
        if sum(r['n'] for r in records) != int(mask.sum()):
            raise RuntimeError('Binning coverage mismatch')
        rows[axis] = records
    ev.dump(DEST/'results.json', dict(selected=len(raw), reference_n=int(mask.sum()), ring_n=int(ring.sum()),
            ring_charge_n=int((ring&charge).sum()), transport_valid_n=int(valid.sum()), rows=rows,
            retained_raw_beta_above_one=int((mask&(raw>1)).sum()),
            retained_new_beta_above_one=int((mask&(new>1)).sum()),
            C12_C13_delay_difference=ev.stats((delays[1]-delays[0])[mask]),
            holdout='original 98 validation runs; previously inspected at low R',
            propagation='fixed C12 expected material delay in inverse beta; no clipping measured beta',
            overall=ev.compare(raw[mask],new[mask],reference[mask])))
    print('EVALUATED selected',len(raw),'reference',int(mask.sum()), flush=True)
    for axis, records in rows.items():
        for r in records:
            print(axis,r['low'],r['high'],'N',r['n'],'gain',r.get('gain_beta_pct'),flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['extract-train','train','extract-validation','evaluate'])
    parser.add_argument('--workers', type=int, default=6)
    args = parser.parse_args()
    if args.stage.startswith('extract-'):
        extract(args.stage.split('-')[1],args.workers)
    elif args.stage=='train':
        train()
    else:
        evaluate()


if __name__ == '__main__':
    main()
