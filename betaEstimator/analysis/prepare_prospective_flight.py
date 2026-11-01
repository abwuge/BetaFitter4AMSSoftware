#!/usr/bin/env python3
"""Select additional measured flight events before opening any new RICH values."""
import argparse
import hashlib
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import uproot
from research_conditional_mc import ROOT, stamp
from research_flight_conditional import write, digest
from rich_transport import track_at_z

OUT = ROOT / 'results/flight_calibration_research/continuation_20261013/prospective_flight'
SOURCE_LIST = ROOT / 'results/tofcalib_data_audit/newtof_files.txt'
PREVIOUS = ROOT / 'results/tofcalib_validation/preparation.json'
RAW_LIST = Path('/lustre02/user/yachen/public/runlist/TofCalib_new/Data/amsroot_list_pass8')
FIELDS = ['run', 'event', 'itrtrack', 'ibetah', 'tk_q', 'tk_qin', 'tk_chis1', 'tk_rigidity1',
          'tk_pos', 'tk_dir', 'tof_hsumh', 'tof_hsumhu', 'tof_tl', 'tof_tl_NewCalib',
          'tof_leng', 'tof_barid', 'tof_adca', 'tof_edep', 'tof_pos']


def exposure():
    OUT.mkdir(exist_ok=True)
    if (OUT / 'exposure.json').exists():
        raise FileExistsError('Exposure inventory already saved')
    old = json.loads(PREVIOUS.read_text())['records']
    pieces = []; evidence = []
    def one(record):
        path = ROOT / 'results/tofcalib_validation/parts' / ('%05d' % record['source']) / 'original.npz'
        with np.load(path, allow_pickle=False) as f:
            runs = np.unique(f['run'])
        return runs, dict(path=str(path), sha256=digest(path), runs=len(runs))
    with ThreadPoolExecutor(max_workers=4) as pool:
        for i, (runs, note) in enumerate(pool.map(one, old)):
            pieces.append(runs); evidence.append(note)
            if (i + 1) % 512 == 0:
                print('FLIGHT_EXPOSURE_RUNS', i + 1, '/', len(old), flush=True)
    extras = [ROOT / 'results/time_local_beta/training.npz',
              ROOT / 'results/flight_calibration_research/mean_loss_likelihood/measured.npz',
              Path('/tmp/betaEstimator-flight-beta-wide/train.npz'),
              Path('/tmp/betaEstimator-flight-beta-wide/validation.npz'),
              Path('/tmp/betaEstimator-flight-mass-beta/training.npz'),
              Path('/tmp/betaEstimator-rich-transport/data_evaluation_input.npz')]
    for path in extras:
        if not path.exists():
            raise FileNotFoundError('Required documented exposure cache missing ' + str(path))
        with np.load(path, allow_pickle=False) as f:
            if 'run' not in f.files:
                raise ValueError('Exposure cache has no run ' + str(path))
            runs = np.unique(f['run'])
        pieces.append(runs); evidence.append(dict(path=str(path), sha256=digest(path), runs=len(runs)))
    runs = np.unique(np.concatenate(pieces)).astype(np.uint32)
    np.savez_compressed(OUT / 'excluded_runs.npz', run=runs)
    write(OUT / 'exposure.json', dict(complete=True, created_utc=stamp(), excluded_run_count=len(runs),
        excluded_source_files=[r['path'] for r in old], exposure_caches=evidence,
        previous_manifest_sha256=digest(PREVIOUS), output_sha256=digest(OUT / 'excluded_runs.npz'),
        source_sha256=digest(Path(__file__)),
        scope='All original carbon candidates in the 4096 historical NEWTOF files, old wide flight training/evaluation, mass-beta training, prior RICH data input, and current timing calibration runs',
        limitation='Covers enumerated persisted betaEstimator histories, not an assertion about every manual study or the upstream global TOF calibration.',
        new_RICH_values_read=False))
    print('FLIGHT_EXPOSURE_COMPLETE', len(runs), flush=True)


def plan():
    if (OUT / 'measurement_protocol.json').exists():
        raise FileExistsError('Prospective measured selection already declared')
    history = json.loads((OUT / 'exposure.json').read_text())
    if digest(OUT / 'excluded_runs.npz') != history['output_sha256']:
        raise ValueError('Excluded run set changed')
    oldfiles = set(history['excluded_source_files'])
    paths = [p for p in SOURCE_LIST.read_text().splitlines() if p and p not in oldfiles]
    paths.sort(key=lambda p: hashlib.sha256(('prospective_flight:' + p).encode()).digest())
    records = [dict(source=i, path=p) for i, p in enumerate(paths)]
    models = {}
    for path in sorted((ROOT / 'results/flight_calibration_research/continuation_20261013/joint_adc/flight').glob('*/freeze.json')):
        models[str(path)] = digest(path)
    write(OUT / 'measurement_protocol.json', dict(created_utc=stamp(), records=records, maximum_candidate_events=50000,
        read_blocks=128, maximum_files=4096,
        source_order='SHA256(prospective_flight:full_path), excluding all previously prepared NEWTOF source files',
        stop_rule='Stop reading complete blocks once measured central-trajectory candidates reach 50000 or 4096 additional source files. Do not inspect RICH counts or outcome errors.',
        selection='Measured carbon, four native fired/used layers, finite non-sentinel native and NewCalib LT times, positive ADC and Edep, finite Tracker L5-L9 geometry; native Kalman rigidity and quality applied after supplement',
        geometry_preselection='Tracker Hermite position at z=-75.1 cm has max(abs(x),abs(y))<17.5 cm; this enriches the NaF path without reading any RICH branch',
        original_timing_slot='tof_tl_NewCalib.reshape(n,-1,4)[:,0,:], previously used method0 LT baseline',
        sampling='Deterministic SHA256(prospective_event:run:event) to the prespecified maximum, after the measured selection',
        target_charge='Same native InitTrkReFit input rule: tk_qin[0,2] when >8.5, otherwise tk_q[1]',
        measurement_fields=FIELDS, reference_fields_read=[], excluded_runs_sha256=history['output_sha256'],
        exposure_sha256=digest(OUT / 'exposure.json'), source_list_sha256=digest(SOURCE_LIST), raw_list_sha256=digest(RAW_LIST),
        source_sha256=digest(Path(__file__)), interpolation_sha256=digest(ROOT / 'analysis/rich_transport.py'),
        frozen_flight_models=models, model_refit_on_new_events=False, new_RICH_values_read=False))
    print('PROSPECTIVE_FLIGHT_PROTOCOL', len(paths), flush=True)


def read_record(record, excluded):
    with uproot.open(record['path']) as f:
        tree = f['amstreea']; raw = tree.arrays(FIELDS, library='np')
    n = len(raw['run']); t = raw['tof_tl'].astype(float)
    new = raw['tof_tl_NewCalib'].reshape(n, -1, 4)[:, 0, :].astype(float)
    adc = raw['tof_adca'].astype(float); x = raw['tof_leng'].astype(float) / 29.9792458
    q = np.where(raw['tk_qin'][:, 0, 2] < 2.5, raw['tk_q'][:, 1], raw['tk_qin'][:, 0, 2])
    used = (raw['tof_hsumh'] == 4) & (raw['tof_hsumhu'] == 4)
    quality = np.isfinite(q) & (q >= 5.5) & (q < 6.5) & used
    quality &= (raw['itrtrack'] >= 0) & (raw['ibetah'] >= 0)
    for value in [t, new]:
        quality &= np.isfinite(value).all(1) & (value != 0).all(1) & (np.abs(value) < 1000).all(1)
    quality &= np.isfinite(adc).all((1, 2)) & (adc > 0).all((1, 2))
    quality &= np.isfinite(raw['tof_edep']).all(1) & (raw['tof_edep'] > 0).all(1)
    quality &= np.isfinite(x).all(1) & (np.ptp(x, axis=1) > 0)
    quality &= ((raw['tof_barid'] >= 0) & (raw['tof_barid'] < 10)).all(1)
    pos = raw['tk_pos'].astype(float); direction = raw['tk_dir'].astype(float)
    geom = np.isfinite(pos[:, 4:9]).all((1, 2)) & np.isfinite(direction[:, 4:9]).all((1, 2))
    geom &= (np.abs(direction[:, 4:9, 2]) > .05).all(1) & (np.diff(pos[:, 4:9, 2], axis=1) < 0).all(1)
    geom &= (pos[:, 4, 2] > 0) & (pos[:, 5, 2] < 0)
    hit = np.zeros(n, bool); good = np.flatnonzero(geom)
    if len(good):
        xy, metric = track_at_z(pos[good], direction[good], np.full(len(good), -75.1))
        hit[good] = np.isfinite(xy).all(1) & (np.max(abs(xy), axis=1) < 17.5)
    overlap = np.isin(raw['run'], excluded)
    take = quality & geom & hit & ~overlap; ix = np.flatnonzero(take)
    out = {k: raw[k][ix] for k in ['run', 'event', 'itrtrack', 'ibetah', 'tk_q', 'tk_qin', 'tk_chis1', 'tk_rigidity1', 'tk_pos', 'tk_dir', 'tof_hsumh', 'tof_hsumhu']}
    out.update(source=np.full(len(ix), record['source'], np.int32), entry=ix, x=x[ix], t=new[ix], t_native=t[ix],
               adc=adc[ix], edep=raw['tof_edep'][ix], bar=raw['tof_barid'][ix], xy=raw['tof_pos'][ix, :, :2])
    return out, dict(record, raw_events=n, detector_valid=int(quality.sum()), central_candidates=int((quality & geom & hit).sum()),
                    excluded_historical_run_events=int((quality & geom & hit & overlap).sum()), accepted=len(ix))


def prepare():
    p = json.loads((OUT / 'measurement_protocol.json').read_text())
    if (OUT / 'candidates.npz').exists():
        raise FileExistsError('Measured candidates already saved')
    for path, key in [(Path(__file__), 'source_sha256'), (OUT / 'excluded_runs.npz', 'excluded_runs_sha256'), (RAW_LIST, 'raw_list_sha256')]:
        if digest(path) != p[key]:
            raise ValueError('Measured preparation input changed')
    with np.load(OUT / 'excluded_runs.npz') as f:
        excluded = f['run']
    pieces = []; notes = []; accepted = 0
    records = p['records'][:p['maximum_files']]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for start in range(0, len(records), p['read_blocks']):
            block = records[start:start + p['read_blocks']]
            for a, note in pool.map(lambda row: read_record(row, excluded), block):
                pieces.append(a); notes.append(note); accepted += len(a['run'])
            print('PROSPECTIVE_FLIGHT_MEASURED', len(notes), '/', p['maximum_files'], 'central_candidates', accepted, flush=True)
            if accepted >= p['maximum_candidate_events']:
                break
    if not accepted:
        raise ValueError('No new measured candidates')
    data = {k: np.concatenate([a[k] for a in pieces]) for k in pieces[0]}
    ids = (data['run'].astype(np.uint64) << np.uint64(32)) | data['event'].astype(np.uint64)
    if len(np.unique(ids)) != len(ids):
        raise ValueError('Duplicated new run/event')
    rank = np.array([int.from_bytes(hashlib.sha256(('prospective_event:%d:%d' % (r, e)).encode()).digest()[:8], 'big')
                     for r, e in zip(data['run'], data['event'])], np.uint64)
    selected = np.sort(np.argsort(rank)[:p['maximum_candidate_events']]); data = {k: v[selected] for k, v in data.items()}
    if np.intersect1d(data['run'], excluded).size:
        raise ValueError('Historical run leakage')
    by_run = {}
    for line in RAW_LIST.read_text().splitlines():
        if line.strip():
            run = int(Path(line.strip()).name.split('.')[0]); by_run.setdefault(run, []).append(line.strip())
    missing = [int(r) for r in np.unique(data['run']) if int(r) not in by_run]
    if missing:
        raise ValueError('No native raw mapping for selected runs ' + str(missing[:10]))
    np.savez_compressed(OUT / 'candidates.npz', **data)
    parts = OUT / 'native'; parts.mkdir(exist_ok=True); jobs = []
    for source in np.unique(data['source']):
        select = np.flatnonzero(data['source'] == source); directory = parts / ('%05d' % source); directory.mkdir(exist_ok=True)
        q = np.where(data['tk_qin'][select, 0, 2] > 8.5, data['tk_qin'][select, 0, 2], data['tk_q'][select, 1])
        rows = np.column_stack([data['run'][select], data['event'][select], data['itrtrack'][select], data['ibetah'][select], q])
        np.savetxt(directory / 'targets.txt', rows, fmt=['%d', '%d', '%d', '%d', '%.9g'])
        raws = [path for run in np.unique(data['run'][select]) for path in by_run[int(run)]]
        (directory / 'raw_files.txt').write_text('\n'.join(raws) + '\n')
        jobs.append(dict(source=int(source), events=len(select), targets_sha256=digest(directory / 'targets.txt'),
                         raw_files_sha256=digest(directory / 'raw_files.txt')))
    write(OUT / 'preparation.json', dict(complete=True, created_utc=stamp(), events=len(selected), parent_candidates=accepted,
        read_sources=len(notes), runs=len(np.unique(data['run'])), historical_run_overlap=0, records=notes, jobs=jobs,
        protocol_sha256=digest(OUT / 'measurement_protocol.json'), candidates_sha256=digest(OUT / 'candidates.npz'),
        stop_rule_reached=accepted >= p['maximum_candidate_events'], new_RICH_values_read=False))
    print('PROSPECTIVE_FLIGHT_READY', len(selected), len(jobs), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('stage', choices=['exposure', 'plan', 'prepare']); args = parser.parse_args()
    globals()[args.stage]()
