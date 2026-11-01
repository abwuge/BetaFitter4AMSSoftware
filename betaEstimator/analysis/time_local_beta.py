#!/usr/bin/env python3
"""Detector-only adjacent-time fits and exhaustive out-of-period prediction.

Stages are separated so no RICH target is opened until every fit is frozen.
The 15-day grid is UTC calendar anchored, not a claim about exact calibration
release boundaries. A full-cycle guard separates targets and training.
"""
import argparse
import json
import multiprocessing as mp
import os
from pathlib import Path
import time

import numpy as np
import uproot

import flight_beta_wide as w
import train_flight_mass_beta as tr
from prepare_tofcalib_validation import OUT as OLD

OUT = w.rt.ROOT / 'results/time_local_beta'
DAY = 86400
PERIOD = 15 * DAY
GUARD = 15 * DAY
WINDOWS = (180, 360, 720, 1080)
MIN_BIN = 100
TRAINING_CACHE = None
FIELDS = list(dict.fromkeys(tr.TRAIN_BRANCHES + ['itrtrack', 'ibetah', 'irich']))
MISSING = {'tk_rigidityk', 'tk_erigidity1', 'tof_bpatt'}


def load(path, keys=None):
    with np.load(str(path), allow_pickle=False) as f:
        return {k: f[k] for k in (f.files if keys is None else keys)}


def detector_preselection(a):
    """Available detector cuts only; no rigidity proxy or RICH selection."""
    t = np.asarray(a['tof_tl'], float)
    x = np.asarray(a['tof_leng'], float) / w.rt.rc.C_LIGHT
    cp = x - x.mean(axis=1, keepdims=True)
    ct = t - t.mean(axis=1, keepdims=True)
    with np.errstate(divide='ignore', invalid='ignore'):
        u = np.sum(cp * ct, axis=1) / np.sum(cp * cp, axis=1)
    chi = a['tk_chis1'][:, 2, 2, 1, :]
    q = np.where(a['tk_qin'][:, 0, 2] < 2.5, a['tk_q'][:, 1], a['tk_qin'][:, 0, 2])
    ok = np.isfinite(q) & (q >= 5.5) & (q < 6.5)
    ok &= np.isfinite(u) & (u > 0) & np.all(np.isfinite(a['tof_pos']), axis=(1, 2))
    ok &= np.all(np.isfinite(a['tk_pos'][:, 1:8]), axis=(1, 2))
    ok &= (a['tk_pos'][:, 4, 2] > 0) & (a['tk_pos'][:, 5, 2] < 0)
    ok &= np.all((a['tof_barid'] >= 0) & (a['tof_barid'] < 10), axis=1)
    ok &= np.all(np.isfinite(a['tof_edep']) & (a['tof_edep'] > 0), axis=1)
    ok &= np.all(np.isfinite(a['tk_ql'][:, 1:8]) & (a['tk_ql'][:, 1:8] > 0), axis=1)
    ok &= np.all(np.isfinite(chi) & (chi > 0) & (chi <= 10), axis=1)
    return ok


def prepare_part(source):
    src = OLD / 'parts' / ('%05d' % source)
    dest = OUT / 'parts' / ('%05d' % source)
    dest.mkdir(parents=True, exist_ok=True)
    note = dest / 'preparation.json'
    signature = dict(original_sha256=w.ev.digest(src / 'original.npz'),
                     extra_sha256=w.ev.digest(src / 'extra.root') if (src / 'extra.root').exists() else None,
                     script_sha256=w.ev.digest(Path(__file__)))
    if note.exists():
        record = json.loads(note.read_text())
        if any(record['signature'][k] != signature[k] for k in ('original_sha256','extra_sha256')) or record['RICH_selection']:
            raise RuntimeError('Prepared input/selection changed')
        return record
    a = load(src / 'original.npz', [k for k in FIELDS if k not in MISSING])
    a = w.ev.sub(a, np.flatnonzero(detector_preselection(a)))
    a['source'] = np.full(len(a['run']), source, np.int32)
    n = len(a['run'])
    if (src / 'extra.root').exists():
        with uproot.open(str(src / 'extra.root')) as f:
            old = f['beta_extra'].arrays(library='np')
    else:
        old_note = json.loads((src / 'supplement.json').read_text())
        if old_note['selected'] != 0 or old_note['requested'] != 0:
            raise RuntimeError('Missing existing nonempty supplement')
        old = {k: np.empty(0, np.uint32) for k in ('run','event','selected','itrtrack','ibetah')}
        old.update(tk_rigidityk=np.empty((0,4,7),np.float32),
                   tk_erigidity1=np.empty((0,3,7),np.float32),tof_bpatt=np.empty(0,np.int32))
    ids = np.rec.fromarrays([a['run'], a['event']], names='run,event')
    extra_ids = np.rec.fromarrays([old['run'], old['event']], names='run,event')
    order = np.argsort(extra_ids)
    index = np.searchsorted(extra_ids[order], ids)
    safe = np.minimum(index, max(0, len(order)-1))
    known = (index < len(order)) & (extra_ids[order][safe] == ids) if len(order) else np.zeros(n,bool)
    for k in MISSING:
        a[k] = np.zeros((n,) + old[k].shape[1:], dtype=old[k].dtype)
        a[k][known] = old[k][order[safe[known]]]
    valid = np.zeros(n, bool)
    valid[known] = old['selected'][order[safe[known]]] == 1
    for k in ('itrtrack', 'ibetah'):
        if np.any(old[k][order[safe[known & valid]]] != a[k][known & valid]):
            raise RuntimeError('Existing native association mismatch')
    a['native_valid'] = valid
    missing = np.flatnonzero(~known)
    a['missing_native'] = ~known
    # Ring index is copied solely to the producer's association interface.
    # It never changes inclusion, weights, measured training fields or loss.
    with (dest / 'targets.txt').open('w') as f:
        for i in missing:
            charge = a['tk_qin'][i, 0, 2] if a['tk_qin'][i, 0, 2] > 8.5 else a['tk_q'][i, 1]
            f.write('%d %d %d %d %d %.9g\n' % (a['run'][i], a['event'][i],
                    a['itrtrack'][i], a['ibetah'][i], a['irich'][i], charge))
    runs = set(int(x) for x in a['run'][missing])
    files = [p for p in (src / 'raw_files.txt').read_text().splitlines()
             if p.strip() and int(Path(p).name.split('.')[0]) in runs]
    (dest / 'raw_files.txt').write_text('\n'.join(files) + '\n')
    np.savez_compressed(str(dest / 'measured.npz'), **a)
    record = dict(complete=True, source=source, selected=n, missing=len(missing),
                  reused=int(known.sum()), RICH_selection=False, signature=signature)
    w.ev.dump(note, record)
    return record


def prepare(workers, files):
    w.check_freeze()
    OUT.mkdir(exist_ok=True)
    config = OUT / 'configuration.json'
    chosen = np.asarray(json.loads((OLD / 'selection.json').read_text())['file_indices'])
    chosen = chosen[np.unique(np.linspace(0, len(chosen)-1, files).astype(int))]
    settings = dict(file_indices=chosen, target_days=15, guard_days=15,
                    training_days_each_side=list(WINDOWS), min_events_per_rigidity_bin=MIN_BIN,
                    sampling='evenly spaced source indices; detector-only carbon selection',
                    exact_calibration_boundaries_available=False,
                    source_freeze_sha256=w.ev.digest(w.DEST / 'freeze.json'),
                    frozen_prediction_sha256=w.ev.digest(OLD / 'predictions.npz'),
                    script_sha256=w.ev.digest(Path(__file__)))
    previous_indices = []
    if config.exists():
        if (OUT / 'plan.json').exists():
            raise RuntimeError('Training folds already fixed')
        previous = json.loads(config.read_text())
        previous_indices = previous.get('native_existing_indices', previous.get('previous_preparation',previous)['file_indices'])
        settings['previous_preparation'] = previous
    settings['native_existing_indices'] = previous_indices
    settings['window_policy'] = 'independent per rigidity block; smallest listed width satisfying detector-only counts'
    w.ev.dump(config, settings)
    records = []
    with mp.Pool(workers) as pool:
        for i, r in enumerate(pool.imap(prepare_part, chosen), 1):
            records.append(r)
            if i % 32 == 0 or i == len(chosen):
                print('Prepared', i, '/', len(chosen), 'new native targets', sum(x['missing'] for x in records), flush=True)
    (OUT / 'condor').mkdir(exist_ok=True)
    pending_indices = [int(x) for x in chosen if int(x) not in set(previous_indices)]
    prefix = 'native_full' if previous_indices else 'native'
    for worker in range(8):
        (OUT / 'condor' / ('%s_%03d.txt' % (prefix, worker))).write_text(
            '\n'.join('%05d' % x for x in pending_indices[worker::8]) + '\n')
    w.ev.dump(OUT / 'preparation.json', dict(complete=True, records=records,
             n=sum(x['selected'] for x in records), missing=sum(x['missing'] for x in records), RICH_selection=False))


def finalize_part(source):
    dest = OUT / 'parts' / ('%05d' % source)
    if (dest / 'training.npz').exists():
        return json.loads((dest / 'finalization.json').read_text())
    a = load(dest / 'measured.npz')
    rows = np.flatnonzero(a.pop('missing_native'))
    note = json.loads((dest / 'preparation.json').read_text())
    if len(rows):
        expected = w.ev.digest(OLD / 'bin/supplement_light')
        if (dest / '.native_complete').read_text().strip() != expected:
            raise RuntimeError('Native completion/binary mismatch')
        with uproot.open(str(dest / 'extra.root')) as f:
            b = f['beta_extra'].arrays(library='np')
        for k in ('run', 'event'):
            if not np.array_equal(a[k][rows], b[k]):
                raise RuntimeError('Incomplete or reordered native identities')
        if not b['found'].all():
            raise RuntimeError('Native event missing')
        valid = b['selected'] == 1
        for k in ('itrtrack', 'ibetah'):
            if not np.array_equal(a[k][rows][valid], b[k][valid]):
                raise RuntimeError('Association mismatch')
        a['native_valid'][rows] = valid
        for k in MISSING:
            a[k][rows] = b[k]
    a = w.ev.sub(a, np.flatnonzero(a.pop('native_valid')))
    w.set_edges(w.R_EDGES)
    a = w.ev.sub(a, np.flatnonzero(w.selection(a)))
    d = w.rt.center_features(a)
    d.update(run=a['run'], event=a['event'], source=a['source'],
             relR=np.abs(a['tk_erigidity1'][:, 2, 1] * d['R']))
    if set(d) != tr.TRAIN_KEYS:
        raise RuntimeError('Answer-bearing training cache')
    np.savez_compressed(str(dest / 'training.npz'), **d)
    r = dict(complete=True, source=source, n=len(d['u']),
             bins=np.bincount(d['bin'], minlength=12), RICH_selection=False,
             cache_sha256=w.ev.digest(dest / 'training.npz'), preparation=note)
    w.ev.dump(dest / 'finalization.json', r)
    return r


def training_mask(run, low, high, days):
    width = days * DAY
    return ((run >= low-GUARD-width) & (run < low-GUARD)) | ((run >= high+GUARD) & (run < high+GUARD+width))


def plan(workers):
    settings = json.loads((OUT / 'configuration.json').read_text())
    pieces, records = [], []
    with mp.Pool(workers) as pool:
        for r in pool.imap(finalize_part, settings['file_indices']):
            records.append(r)
            pieces.append(load(OUT / 'parts' / ('%05d' % r['source']) / 'training.npz'))
    a = w.ev.concatenate(pieces)
    ids = np.rec.fromarrays([a['run'], a['event']], names='run,event')
    if len(np.unique(ids)) != len(ids):
        raise RuntimeError('Training duplicates')
    np.savez_compressed(str(OUT / 'training.npz'), **a)
    # Only identity and time, not RICH values or residuals, are read for folds.
    run = load(OLD / 'predictions.npz', ['run'])['run'].astype(np.int64)
    anchor = int(run.min() // DAY * DAY)
    fold_ids = (run-anchor) // PERIOD
    folds = []
    for f in np.unique(fold_ids):
        low, high = anchor+int(f)*PERIOD, anchor+(int(f)+1)*PERIOD
        blocks = []
        for i, edges in enumerate(w.BLOCKS):
            choices = []
            for days in WINDOWS:
                train = training_mask(a['run'], low, high, days)
                counts = np.histogram(a['R'][train], bins=edges)[0]
                choices.append(dict(days=days, n=int(counts.sum()), bins=counts))
                if counts.min() >= MIN_BIN:
                    break
            choice = choices[-1]
            if min(choice['bins']) < MIN_BIN:
                raise RuntimeError('Insufficient measured training coverage for fold %d block %d: %s' % (f, i, choice))
            blocks.append(dict(block=i, training_days=choice['days'], train_n=choice['n'], train_bins=choice['bins'], choices=choices))
        max_days = max(b['training_days'] for b in blocks)
        target_runs = np.unique(run[fold_ids == f])
        train_runs = np.unique(a['run'][training_mask(a['run'], low, high, max_days)])
        if np.intersect1d(target_runs, train_runs).size:
            raise RuntimeError('Target run in training')
        folds.append(dict(fold=int(f), start=low, stop=high, blocks=blocks,
                          target_n=int((fold_ids == f).sum()), train_run_n=len(train_runs), target_run_n=len(target_runs)))
    w.ev.dump(OUT / 'plan.json', dict(complete=True, anchor=anchor, folds=folds,
             target_n=len(run), train_n=len(a['run']), guard_days=15,
             training_fields=sorted(a), RICH_selection=False, RICH_used_for_fitting=False,
             training_sha256=w.ev.digest(OUT / 'training.npz'),
             configuration_sha256=w.ev.digest(OUT / 'configuration.json'),
             dependencies={str(Path(p).resolve()): w.ev.digest(p) for p in
                           [__file__, tr.__file__, w.__file__, w.rt.__file__, w.rt.rc.__file__, w.rt.er.__file__]}))
    (OUT / 'folds').mkdir(exist_ok=True)
    batch_jobs = int(settings.get('batch_jobs', 50))
    for worker in range(batch_jobs):
        (OUT / 'condor' / ('fit_%03d.txt' % worker)).write_text('\n'.join(str(x['fold']) for x in folds[worker::batch_jobs])+'\n')
    print('PLANNED', len(folds), 'folds;', len(a['run']), 'detector-only training events', flush=True)


def fit_one(fold_id):
    global TRAINING_CACHE
    plan_data = json.loads((OUT / 'plan.json').read_text())
    for p, digest in plan_data['dependencies'].items():
        if w.ev.digest(p) != digest:
            raise RuntimeError('Training dependency changed')
    record = next(f for f in plan_data['folds'] if f['fold'] == fold_id)
    dest = OUT / 'folds' / ('%04d' % fold_id)
    dest.mkdir(exist_ok=True)
    if (dest / 'freeze.json').exists():
        return
    if TRAINING_CACHE is None:
        a = load(OUT / 'training.npz')
        if set(a) != tr.TRAIN_KEYS:
            raise RuntimeError('Unexpected training fields')
        if w.ev.digest(OUT / 'training.npz') != plan_data['training_sha256']:
            raise RuntimeError('Training cache changed')
        TRAINING_CACHE = a
    a = TRAINING_CACHE
    max_days = max(b['training_days'] for b in record['blocks'])
    np.savez_compressed(str(dest / 'training_runs.npz'), run=np.unique(a['run'][training_mask(a['run'], record['start'], record['stop'], max_days)]))
    blocks = []
    for i, edges in enumerate(w.BLOCKS):
        w.set_edges(edges)
        days = record['blocks'][i]['training_days']
        selected = training_mask(a['run'], record['start'], record['stop'], days)
        d = w.ev.sub(a, selected & (a['R'] >= edges[0]) & (a['R'] < edges[-1]))
        d['bin'] = np.searchsorted(edges, d['R'], side='right')-1
        counts = np.bincount(d['bin'], minlength=4)
        if counts.min() < MIN_BIN:
            raise RuntimeError('Training coverage differs from fixed plan')
        raw = tr.fit_raw(d)
        state, matrix = tr.features(d)
        try:
            params, fitted = tr.fit_correction(d, matrix, raw, raw['fraction'])
        except RuntimeError as error:
            if 'ITERATION' not in str(error).upper():
                raise
            # Solver budget only: unchanged objective, bounds, features and tolerance.
            original_minimize = tr.minimize
            def extended_minimize(*args, **kwargs):
                options = dict(kwargs.get('options', {}))
                options['maxiter'] = 20000
                kwargs['options'] = options
                return original_minimize(*args, **kwargs)
            try:
                tr.minimize = extended_minimize
                params, fitted = tr.fit_correction(d, matrix, raw, raw['fraction'])
                fitted['solver_retry_maxiter'] = 20000
            finally:
                tr.minimize = original_minimize
        if not np.isfinite(params).all() or not np.isfinite(fitted['objective']):
            raise RuntimeError('Non-finite local fitted parameters/objective')
        path = dest / ('block_%d_model.json' % i)
        w.ev.dump(path, dict(state, parameters=params, rigidity_edges_GV=edges))
        replay = tr.feature_matrix(d, w.rt.arrays(json.loads(path.read_text())))
        error = float(np.max(np.abs(matrix.dot(params[:matrix.shape[1]])-replay.dot(params[:replay.shape[1]]))))
        if error > 1e-12:
            raise RuntimeError('Serialization differs')
        blocks.append(dict(fitted, n=len(d['u']), counts=counts, edges=edges,
                           raw_spectrum=raw, training_days=days, serialization_max_difference=error))
        print('FIT', fold_id, 'block', i, 'N', len(d['u']), 'iterations', fitted['iterations'], flush=True)
    w.ev.dump(dest / 'training.json', dict(complete=True, blocks=blocks, fold=record,
             RICH_used=False, MC_used=False, RICH_selection=False,
             runtime=dict(host=os.uname().nodename, numpy=np.__version__, scipy=__import__('scipy').__version__)))
    w.ev.dump(dest / 'freeze.json', dict(utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
             plan_sha256=w.ev.digest(OUT / 'plan.json'), training_sha256=plan_data['training_sha256'],
             model_sha256={str(i): w.ev.digest(dest / ('block_%d_model.json' % i)) for i in range(3)},
             training_runs_sha256=w.ev.digest(dest / 'training_runs.npz'), fold=record))
    if (dest / 'failure.json').exists():
        (dest / 'failure.json').unlink()



def fit_safely(fold_id):
    import traceback
    try:
        fit_one(fold_id)
        return True
    except Exception:
        dest = OUT / 'folds' / ('%04d' % fold_id)
        dest.mkdir(exist_ok=True)
        w.ev.dump(dest / 'failure.json', dict(complete=False, fold=fold_id,
                  traceback=traceback.format_exc(), RICH_opened=False))
        return False

def freeze_all():
    plan_data = json.loads((OUT / 'plan.json').read_text())
    snapshots = {}
    for f in plan_data['folds']:
        dest = OUT / 'folds' / ('%04d' % f['fold'])
        freeze = json.loads((dest / 'freeze.json').read_text())
        if freeze['plan_sha256'] != w.ev.digest(OUT / 'plan.json'):
            raise RuntimeError('Fit plan changed')
        for i, digest in freeze['model_sha256'].items():
            if w.ev.digest(dest / ('block_%s_model.json' % i)) != digest:
                raise RuntimeError('Local model changed')
        snapshots[str(f['fold'])] = w.ev.digest(dest / 'freeze.json')
    w.ev.dump(OUT / 'freeze.json', dict(complete=True, folds=snapshots,
             plan_sha256=w.ev.digest(OUT / 'plan.json'), RICH_used_for_fitting=False))
    return plan_data


def evaluate():
    # RICH becomes accessible only after ALL candidate models have been frozen.
    plan_data = freeze_all()
    w.check_freeze()
    p = load(OLD / 'predictions.npz')
    c = load(OLD / 'diagnostic_covariates.npz', ['run', 'event', 'source', 'continuous', 'bars', 'relR'])
    for k in ('run', 'event', 'source'):
        if not np.array_equal(p[k], c[k]):
            raise RuntimeError('Validation features reordered')
    fold_ids = ((p['run'].astype(np.int64)-plan_data['anchor']) // PERIOD).astype(np.int32)
    corrected = np.full(len(p['run']), np.nan)
    target_count = np.zeros(len(p['run']), np.int8)
    periods = []
    for f in plan_data['folds']:
        rows = np.flatnonzero(fold_ids == f['fold'])
        dest = OUT / 'folds' / ('%04d' % f['fold'])
        train_runs = load(dest / 'training_runs.npz')['run']
        if np.intersect1d(train_runs, p['run'][rows]).size:
            raise RuntimeError('Validation run leakage')
        for i, edges in enumerate(w.BLOCKS):
            w.set_edges(edges)
            local_rows = rows[(p['R'][rows] >= edges[0]) & (p['R'][rows] < edges[-1])]
            d = dict(u=p['raw_u'][local_rows], R=p['R'][local_rows],
                     continuous=c['continuous'][local_rows].astype(float), bars=c['bars'][local_rows])
            d['bin'] = np.searchsorted(edges, d['R'], side='right')-1
            model = w.rt.arrays(json.loads((dest / ('block_%d_model.json' % i)).read_text()))
            x = tr.feature_matrix(d, model)
            corrected[local_rows] = d['u']-x.dot(model['parameters'][:x.shape[1]])
        target_count[rows] += 1
    if not np.isfinite(corrected).all() or not np.all(target_count == 1):
        raise RuntimeError('Incomplete or multiply predicted target events')
    new = 1/(corrected+p['delay_u_C12'])
    mask = p['reference_mask'] & np.isfinite(new) & (new > 0)
    if not np.array_equal(mask, p['reference_mask']):
        raise RuntimeError('Invalid local prediction; cannot silently change paired evaluation sample')
    for f in plan_data['folds']:
        rows = mask & (fold_ids == f['fold'])
        r = dict(f, reference_n=int(rows.sum()))
        if rows.sum() >= 2:
            r['local'] = w.ev.compare(p['raw_beta_rich'][rows], new[rows], p['rich_beta'][rows])
            r['frozen'] = w.ev.compare(p['raw_beta_rich'][rows], p['new_beta_rich'][rows], p['rich_beta'][rows])
        periods.append(r)
    years = p['run'].astype('datetime64[s]').astype('datetime64[Y]').astype(int)+1970
    yearly = []
    for y in np.unique(years):
        rows = mask & (years == y)
        r = dict(year=int(y), n=int(rows.sum()),
                 local=w.ev.compare(p['raw_beta_rich'][rows], new[rows], p['rich_beta'][rows]),
                 frozen=w.ev.compare(p['raw_beta_rich'][rows], p['new_beta_rich'][rows], p['rich_beta'][rows]))
        r['bootstrap'] = w.ev.bootstrap(p['raw_beta_rich'][rows], new[rows], p['rich_beta'][rows], p['run'][rows])
        yearly.append(r)
        print('YEAR', y, 'N', r['n'], 'local', r['local']['gain_beta_pct'], 'frozen', r['frozen']['gain_beta_pct'], flush=True)
    from evaluate_tofcalib_validation import R_EDGES, BETA_EDGES
    axes = {}
    for axis, values, edges in [('rigidity', p['R'], R_EDGES), ('beta', p['rich_beta'], BETA_EDGES)]:
        records = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            rows = mask & (values >= lo) & (values < hi)
            r = dict(low=lo, high=hi, n=int(rows.sum()))
            if rows.sum() >= 2:
                r.update(w.ev.compare(p['raw_beta_rich'][rows], new[rows], p['rich_beta'][rows]))
                r['frozen'] = w.ev.compare(p['raw_beta_rich'][rows], p['new_beta_rich'][rows], p['rich_beta'][rows])
                r['x'] = float(np.median(values[rows]))
                r['bootstrap'] = w.ev.bootstrap(p['raw_beta_rich'][rows], new[rows], p['rich_beta'][rows], p['run'][rows])
            records.append(r)
        if sum(r['n'] for r in records) != int(mask.sum()):
            raise RuntimeError('Axis coverage differs')
        axes[axis] = records
    np.savez_compressed(str(OUT / 'predictions.npz'), **dict(p, local_corrected_u=corrected, local_beta_rich=new, fold=fold_ids))
    w.ev.dump(OUT / 'results.json', dict(complete=True, reference_n=int(mask.sum()), period_rows=periods,
             year_rows=yearly, rows=axes, RICH_used_for_fitting=False, RICH_selection_for_training=False,
             target_run_leakage=0, each_event_predicted_once=True, guard_days=15,
             overall=w.ev.compare(p['raw_beta_rich'][mask], new[mask], p['rich_beta'][mask]),
             frozen_overall=w.ev.compare(p['raw_beta_rich'][mask], p['new_beta_rich'][mask], p['rich_beta'][mask]),
             prediction_sha256=w.ev.digest(OUT / 'predictions.npz'),
             limitations=['15-day UTC grid does not establish exact native calibration boundaries',
                          'conditional run-bootstrap intervals do not refit neighbouring models',
                          'validation continuous features come from existing float32 diagnostic cache',
                          'window increases depend on training counts only; edge folds can be one-sided',
                          'existing reference has previously been inspected; this is a follow-up comparison, not a pristine test']))
    plot(periods, yearly, axes)


def plot(periods, yearly, axes):
    from plot_time_local_beta import plot_trends
    plot_trends(OUT, dict(period_rows=periods, year_rows=yearly, rows=axes))


def main():
    global TRAINING_CACHE
    ap = argparse.ArgumentParser()
    ap.add_argument('stage', choices=['prepare','plan','fit','evaluate'])
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--files', type=int, default=4096)
    ap.add_argument('--fold', type=int)
    ap.add_argument('--batch', type=Path)
    args = ap.parse_args()
    if args.stage == 'prepare': prepare(args.workers, args.files)
    elif args.stage == 'plan': plan(args.workers)
    elif args.stage == 'fit':
        if args.batch:
            tasks = [int(x) for x in args.batch.read_text().splitlines() if x.strip()]
            plan_data = json.loads((OUT / 'plan.json').read_text())
            TRAINING_CACHE = load(OUT / 'training.npz')
            if set(TRAINING_CACHE) != tr.TRAIN_KEYS or w.ev.digest(OUT / 'training.npz') != plan_data['training_sha256']:
                raise RuntimeError('Training cache schema/hash changed')
            with mp.Pool(args.workers) as pool:
                failed = 0
                for i, success in enumerate(pool.imap_unordered(fit_safely, tasks), 1):
                    failed += int(not success)
                    print('FINISHED', i, '/', len(tasks), 'failed', failed, flush=True)
                if failed:
                    raise RuntimeError('%d folds failed; see per-fold failure.json' % failed)
        else:
            fit_one(args.fold)
    else: evaluate()


if __name__ == '__main__':
    main()
