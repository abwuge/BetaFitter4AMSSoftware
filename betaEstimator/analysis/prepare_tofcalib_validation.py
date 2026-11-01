#!/usr/bin/env python3
"""Prepare a time-spanning NEWTOF holdout and exact native-event targets."""
import argparse
import json
import multiprocessing as mp
from pathlib import Path

import numpy as np
import uproot
import flight_beta_wide as w

OUT = w.rt.ROOT / 'results/tofcalib_validation'
EXTRA = {'tk_rigidityk', 'tk_erigidity1', 'tof_bpatt', 'rich_itrtrack', 'rich_ebeta1'}
FIELDS = list(dict.fromkeys(w.ev.MEASURED + w.ev.SELECTION + ['rich_q', 'ibetah']))
RAW_LIST = Path('/lustre02/user/yachen/public/runlist/TofCalib_new/Data/amsroot_list_pass8')
TRAIN_RUNS = None
RAW_BY_RUN = None


def initialize():
    global TRAIN_RUNS, RAW_BY_RUN
    with np.load(str(w.CACHE / 'train.npz'), allow_pickle=False) as f:
        TRAIN_RUNS = np.unique(f['run'])
    RAW_BY_RUN = {}
    for line in RAW_LIST.read_text().splitlines():
        path = line.strip()
        if path:
            run = int(Path(path).name.split('.')[0])
            RAW_BY_RUN.setdefault(run, []).append(path)


def prepare(task):
    source, path = task
    dest = OUT / 'parts' / ('%05d' % source)
    dest.mkdir(parents=True, exist_ok=True)
    note = dest / 'preparation.json'
    if note.exists():
        old = json.loads(note.read_text())
        if old.get('complete') and old.get('path') == path:
            return old
    with uproot.open(path) as f:
        tree = f['amstreea']
        a = tree.arrays(['run', 'event', 'tk_qin', 'tk_q'], library='np')
        q = np.where(a['tk_qin'][:, 0, 2] < 2.5, a['tk_q'][:, 1], a['tk_qin'][:, 0, 2])
        mask = np.isfinite(q) & (q >= 5.5) & (q < 6.5) & ~np.isin(a['run'], TRAIN_RUNS)
        rows = np.flatnonzero(mask)
        # No beta or RICH performance cut is used to prepare targets. Exact
        # frozen selections are applied only after supplementing measurements.
        fields = [k for k in FIELDS if k not in EXTRA]
        original = tree.arrays(fields, library='np')
        data = w.ev.sub(original, rows)
        data['source'] = np.full(len(rows), source, np.int32)
        data['original_entry'] = rows
        raw = int(tree.num_entries)
    runs = np.unique(data['run'])
    missing = [int(r) for r in runs if int(r) not in RAW_BY_RUN]
    if missing:
        raise RuntimeError('Raw AMS mapping missing runs: %s' % missing[:10])
    native_inputs = [p for r in runs for p in RAW_BY_RUN[int(r)]]
    (dest / 'raw_files.txt').write_text('\n'.join(native_inputs) + '\n')
    with (dest / 'targets.txt').open('w') as f:
        for run, event in zip(data['run'], data['event']):
            f.write('%d %d\n' % (run, event))
    np.savez_compressed(str(dest / 'original.npz'), **data)
    record = dict(complete=True, path=path, source=source, raw=raw,
                  candidates=len(rows), actual_runs=len(runs), raw_files=len(native_inputs),
                  excluded_training_runs=True, beta_or_rich_target_selection=False)
    w.ev.dump(note, record)
    return record


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--files', type=int, default=4096)
    p.add_argument('--workers', type=int, default=8)
    p.add_argument('--limit', type=int)
    args = p.parse_args()
    w.check_freeze()
    OUT.mkdir(exist_ok=True)
    all_files = (w.rt.ROOT / 'results/tofcalib_data_audit/newtof_files.txt').read_text().splitlines()
    # Even coverage across file start times, chosen before looking at outcomes.
    all_files.sort(key=lambda x: (int(Path(x).name.split('_')[0]), x))
    chosen = np.unique(np.linspace(0, len(all_files)-1, min(args.files, len(all_files))).astype(int))
    tasks = [(int(i), all_files[int(i)]) for i in chosen]
    w.ev.dump(OUT / 'selection.json', dict(total_source_files=len(all_files),
              requested_source_files=len(tasks), file_indices=chosen,
              sampling='evenly spaced file start times; no outcome-based selection',
              exclude_all_actual_training_runs=True, frozen_model_sha256=w.ev.digest(w.DEST / 'freeze.json')))
    if args.limit:
        tasks = tasks[:args.limit]
    records = []
    with mp.Pool(args.workers, initializer=initialize) as pool:
        for i, record in enumerate(pool.imap_unordered(prepare, tasks), 1):
            records.append(record)
            if i % 10 == 0 or i == len(tasks):
                print('Prepared', i, '/', len(tasks), 'candidates', sum(r['candidates'] for r in records), flush=True)
    w.ev.dump(OUT / 'preparation.json', dict(complete=len(tasks)==len(chosen),
              completed_files=len(records), candidates=sum(r['candidates'] for r in records), records=records))


if __name__ == '__main__':
    main()
