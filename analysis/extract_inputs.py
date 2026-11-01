#!/usr/bin/env python3
"""Extract detector-only arrays. Never requests simulation branches."""
import argparse
import hashlib
import json
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import uproot

CHARGES = (1, 2, 6, 8)
BRANCHES = ['run', 'event', 'tof_betah', 'tof_tl', 'tof_leng', 'tof_barid',
            'tof_pos', 'tof_dir', 'tof_edep', 'tof_bpatt', 'tk_edep',
            'tk_pos', 'tk_q', 'tk_qin']
ARRAY_KEYS = ('run', 'event', 'tof_betah', 'tof_tl', 'tof_leng', 'tof_barid',
              'tof_edep', 'tof_bpatt', 'tk_edep', 'position', 'direction',
              'charge', 'source_id', 'source_entry')


def parse_sources(path):
    records, seen = [], set()
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        role, z, name = line.split()
        if role not in ('train', 'validation') or name in seen:
            raise ValueError('Invalid role or duplicate source: ' + line)
        seen.add(name)
        records.append((role, int(z), name))
    return records


def detector_selection(a):
    q = np.where(a['tk_qin'][:, 0, 2] < 2.5, a['tk_q'][:, 1], a['tk_qin'][:, 0, 2])
    charge = np.where(np.isfinite(q), q + .5, -1).astype(np.int32)
    keep = np.isin(charge, CHARGES)
    keep &= np.all(np.isfinite(a['tk_edep'][:, 1:8]) & (a['tk_edep'][:, 1:8] > 0), axis=1)
    keep &= np.all(np.isfinite(a['tk_pos'][:, 1:8]), axis=(1, 2))
    keep &= (a['tk_pos'][:, 4, 2] > 0) & (a['tk_pos'][:, 5, 2] < 0)
    return charge, keep


def load_raw(path):
    with uproot.open(path) as f:
        a = f['amstreea'].arrays(BRANCHES, library='np')
    charge, keep = detector_selection(a)
    out = {k: a[k][keep] for k in ARRAY_KEYS if k in a}
    out['position'] = a['tof_pos'][keep, :, :2]
    out['direction'] = a['tof_dir'][keep, 0]
    out['charge'] = charge[keep]
    out['source_entry'] = np.flatnonzero(keep)
    return out, len(charge)


def load_arrays(path):
    with np.load(str(path), allow_pickle=False) as f:
        if set(f.files) != set(ARRAY_KEYS):
            raise ValueError('Expected detector-only cache; regenerate ' + str(path))
        return {k: f[k] for k in ARRAY_KEYS}


def read_record(record):
    source_id, (role, source_z, path) = record
    d, total = load_raw(path)
    return source_id, role, path, d, total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sources', type=Path, default=Path('config/data_sources.conf'))
    ap.add_argument('--roles', nargs='+', default=['train', 'validation'], choices=['train', 'validation'])
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    stores = defaultdict(lambda: defaultdict(list))
    audit, run_roles, event_ids = [], {}, set()
    records = [(i, record) for i, record in enumerate(parse_sources(args.sources), 1)
               if record[0] in args.roles]
    executor = ProcessPoolExecutor(max_workers=4)
    for source_id, role, path, d, total in executor.map(read_record, records):
        for run in np.unique(d['run']):
            if run_roles.setdefault(int(run), role) != role:
                raise ValueError('Run occurs in training and validation: %s' % run)
        ids = list(zip(d['run'].tolist(), d['event'].tolist()))
        if len(set(ids)) != len(ids) or any(key in event_ids for key in ids):
            raise ValueError('Duplicate (run,event) in ' + path)
        event_ids.update(ids)
        d['source_id'] = np.full(len(d['charge']), source_id, dtype=np.int32)
        for z in CHARGES:
            mask = d['charge'] == z
            for k, v in d.items():
                stores[(role, z)][k].append(v[mask])
        audit.append(dict(role=role, source_id=source_id, path=path,
                          raw_entries=total, selected_entries=len(d['charge']),
                          bytes=Path(path).stat().st_size))
        if len(audit) % 20 == 0:
            print('read %d sources' % len(audit), flush=True)
    executor.shutdown()
    counts = {}
    for (role, z), store in sorted(stores.items()):
        out = {k: np.concatenate(v) for k, v in store.items()}
        dest = args.output / ('%s_Z%d.npz' % (role, z))
        np.savez_compressed(str(dest), **out)
        counts[dest.name] = len(out['charge'])
        print('%s rows=%d' % (dest.name, len(out['charge'])), flush=True)
    (args.output / 'extract_audit.json').write_text(json.dumps(dict(
        sources_sha256=hashlib.sha256(args.sources.read_bytes()).hexdigest(),
        branches=BRANCHES, counts=counts, sources=audit,
        disjoint_runs=True, unique_events=True), indent=2))


if __name__ == '__main__':
    main()
