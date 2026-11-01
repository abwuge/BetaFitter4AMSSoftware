#!/usr/bin/env python3
"""Check completed artifacts and add paired, fixed-model mass-separation CIs."""
import csv
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np
import uproot

root = Path(__file__).resolve().parent
def read(name):
    return list(csv.DictReader((root/name).open(), delimiter='\t'))
def width(v):
    q = np.percentile(v, [16, 84]); return float((q[1]-q[0])/2)
def separation(v, tag):
    a, b = v[tag == 0], v[tag == 1]
    return float((np.median(b)-np.median(a))/np.sqrt(width(a)**2+width(b)**2))
with uproot.open(str(root/'predictions.root')) as f:
    a = f['predictions'].arrays(['Z','run','event','tag','rawU','calibratedU','referenceU',
        'mass2Raw','mass2Calibrated','mass2Reference','R','Rinjected',
        'mass2RigidityInjectionBefore','mass2RigidityInjectionAfter'], library='np')
frozen = json.loads((root/'frozen_models.json').read_text())
if frozen['script_sha256'] != hashlib.sha256((root/'fit_mass_study.py').read_bytes()).hexdigest():
    raise ValueError('Frozen training script hash mismatch')
ids = np.rec.fromarrays([a['run'], a['event']])
if len(ids) != len(np.unique(ids)):
    raise ValueError('Duplicate predictions')
for k in ['rawU','calibratedU','referenceU','mass2Raw','mass2Calibrated','mass2Reference']:
    if not np.isfinite(a[k]).all():
        raise ValueError('Nonfinite reconstruction '+k)
sources = read('flight.root.sources.tsv')
report = {'truth_access': False, 'prediction_entries': len(ids),
          'source_count': len(sources), 'raw_source_entries': sum(int(r['raw']) for r in sources),
          'extracted_entries': sum(int(r['kept']) for r in sources),
          'training_code_hash_matches': True, 'predictions_unique': True,
          'gradient_checks': frozen['gradient_checks'], 'rich_training_counts': frozen['rich_training_counts'],
          'charges': {}}
boot = read('bootstrap.tsv')
for z in [1,2]:
    q = a['Z'] == z
    tagq = q & (a['tag'] >= 0)
    allruns = np.unique(a['run'][q])
    rows = {int(r): np.flatnonzero(tagq & (a['run'] == r)) for r in allruns}
    rng = np.random.RandomState(14031+z)
    differences = []
    for i in range(200):
        draws = rng.choice(allruns, len(allruns), replace=True)
        ix = np.concatenate([rows[int(r)] for r in draws])
        labels = a['tag'][ix]
        if np.sum(labels == 0) < 10 or np.sum(labels == 1) < 10:
            continue
        before = separation(a['mass2Raw'][ix], labels)
        after = separation(a['mass2Calibrated'][ix], labels)
        differences.append(after/before-1)
    bw = np.array([float(r['rich_width_narrowing']) for r in boot if int(r['Z']) == z])
    bu = np.array([float(r['auc_difference']) for r in boot if int(r['Z']) == z])
    p = np.flatnonzero(q)[np.argmax(a['rawU'][q])]
    report['charges'][str(z)] = {
        'n_holdout': int(q.sum()), 'n_runs': len(allruns),
        'rich_width_narrowing_ci95': np.percentile(bw[np.isfinite(bw)], [2.5,97.5]).tolist(),
        'aggregate_tag_auc_difference_ci95': np.percentile(bu[np.isfinite(bu)], [2.5,97.5]).tolist(),
        'tag_mass2_separation_relative_change_ci95': np.percentile(differences,[2.5,97.5]).tolist(),
        'nonpositive_raw_u': int(np.sum(a['rawU'][q] <= 0)),
        'nonpositive_calibrated_u': int(np.sum(a['calibratedU'][q] <= 0)),
        'raw_u_above_2': int(np.sum(a['rawU'][q] > 2)),
        'raw_u_above_5': int(np.sum(a['rawU'][q] > 5)),
        'raw_u_percentiles': np.percentile(a['rawU'][q],[.1,1,50,99,99.9,100]).tolist(),
        'correction_percentiles': np.percentile((a['rawU']-a['calibratedU'])[q],[0,1,50,99,100]).tolist(),
        'largest_u_event': {'run': int(a['run'][p]), 'event': int(a['event'][p]),
                            'R':float(a['R'][p]), 'rawU':float(a['rawU'][p])}}
info = subprocess.check_output(['pdfinfo',str(root/'mass-study.pdf')]).decode('utf8')
report['pdf_pages'] = int(next(x.split(':',1)[1] for x in info.splitlines() if x.startswith('Pages:')))
if report['pdf_pages'] != 14:
    raise ValueError('Unexpected PDF page count')
gaussian = read('gaussian_core.tsv')
report['gaussian_core_fits'] = len(gaussian)
if len(gaussian) != 8 or any(int(x['fit_status']) != 0 for x in gaussian):
    raise ValueError('Gaussian core diagnostic fit failed')
(root/'qa.json').write_text(json.dumps(report,indent=2,allow_nan=False))
print(json.dumps(report,indent=2,allow_nan=False))
