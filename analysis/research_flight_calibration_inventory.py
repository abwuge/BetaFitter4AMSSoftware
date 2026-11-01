#!/usr/bin/env python3
"""Inventory measured inputs for a flight-only calibration study."""
import json
import hashlib
from pathlib import Path
import numpy as np
import uproot

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/flight_calibration_research'
MEASURED = ['run', 'event', 'version', 'tof_tl', 'tof_leng', 'tof_etl',
            'tof_pos', 'tof_barid', 'tof_edep', 'tof_rtime', 'tof_Aadc', 'tof_adca',
            'tof_rtime_ht', 'tof_rtime_sht', 'tof_sbit', 'tk_edep',
            'tk_ql', 'tk_qlqs', 'tk_pos', 'tk_dir']


def describe(values):
    a = np.asarray(values)
    r = dict(shape=list(a.shape), dtype=str(a.dtype), finite_fraction=float(np.isfinite(a).mean()))
    if a.size:
        v = a[np.isfinite(a)]
        r['quantiles'] = np.percentile(v, [0, 5, 50, 95, 100]).tolist() if len(v) else []
        r['nonzero_fraction'] = float((a != 0).mean())
        r['minus_one_fraction'] = float((a == -1).mean())
    return r


def main():
    OUT.mkdir(exist_ok=True)
    plan = json.loads((ROOT/'results/time_local_beta/configuration.json').read_text())
    indices = plan['file_indices']
    selected = [int(indices[j]) for j in np.unique(np.linspace(0, len(indices)-1, 9).astype(int))]
    records = json.loads((ROOT/'results/tofcalib_validation/preparation.json').read_text())['records']
    paths = {int(r['source']): r['path'] for r in records}
    result = dict(RICH_values_read=False, MC_values_read=False, sample_policy='nine evenly spaced prepared source indices; at most first 5000 original entries', records=[])
    for source in selected:
        item = dict(source=source, path=paths[source])
        for role, path in [('training_measured', ROOT/'results/time_local_beta/parts'/('%05d'%source)/'measured.npz'),
                           ('original_cache', ROOT/'results/tofcalib_validation/parts'/('%05d'%source)/'original.npz')]:
            with np.load(str(path), allow_pickle=False) as f:
                item[role] = dict(available_measured=[k for k in MEASURED if k in f.files], missing_measured=[k for k in MEASURED if k not in f.files])
        with uproot.open(paths[source]) as f:
            tree = f['amstreea']
            keys = set(tree.keys())
            item['other_response_branch_names'] = sorted(k for k in keys if any(s in k.lower() for s in ['adc', 'edep', 'rtime', 'qlpm', 'dedx', 'amplitude']))
            chosen = [k for k in MEASURED if k in keys]
            a = tree.arrays(chosen, entry_stop=min(5000, tree.num_entries), library='np')
            item['raw_entries'] = int(tree.num_entries)
            item['raw_missing'] = [k for k in MEASURED if k not in keys]
            item['raw_columns'] = {k: describe(a[k]) for k in chosen}
        result['records'].append(item)
        print('source', source, 'raw missing', item['raw_missing'], 'response branches', item['other_response_branch_names'], flush=True)
        print('rtime', {k:item['raw_columns'][k] for k in ('tof_rtime','tof_rtime_ht','tof_rtime_sht','tof_adca') if k in item['raw_columns']}, flush=True)
    result['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    path = OUT/'input_inventory.json'
    path.write_text(json.dumps(result, indent=2)+'\n')
    print('COMPLETE', path, flush=True)

if __name__ == '__main__':
    main()
