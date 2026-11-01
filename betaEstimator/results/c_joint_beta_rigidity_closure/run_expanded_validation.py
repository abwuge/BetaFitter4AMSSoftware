#!/usr/bin/env python3
"""Evaluate frozen models on additional whole files, retaining all selected events."""
import csv
import glob
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / 'c_isotope_proportion_closure'
OUT = HERE / 'expanded_validation'
sys.path.insert(0, str(BASE))
import run_closure as rc
import run_tracker_multifit_probe as rt
import run_joint_identifiability_scan as js

METHODS = ['raw', 'beta_only_frozen', 'joint_unanchored', 'joint_anchored']


def arrays(obj):
    if isinstance(obj, dict):
        return {k: arrays(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return np.asarray(obj)
    return obj


def tracker_features(d, t, state):
    candidates, names = rt.candidate_matrix(t)
    raw, names = rt.raw_quality_features(d, t, candidates, names)
    raw = raw[:, [not n.startswith('geometry_') for n in names]]
    x = rt.apply_transform_state(raw, state['raw_state'])[:, 1:]
    pc = x.dot(state['pca_components'].T) / state['pca_scales']
    return np.column_stack([np.ones(len(pc)), pc])


def anchor(x, d, state):
    x = (x-state['bin_centers'][d['bin']]).dot(state['null_projection']) / state['scale']
    return np.column_stack([np.zeros(len(x)), x])


def predict(d, t, unanchored, anchored):
    xb = rc.apply_feature_state(d, unanchored['beta_feature_state'])
    xr = tracker_features(d, t, unanchored['tracker_feature_state'])
    ab = anchor(rc.apply_feature_state(d, anchored['beta_base_state'])[:, 1:], d, anchored['beta_anchor_state'])
    ar = anchor(tracker_features(d, t, anchored['rigidity_base_state'])[:, 1:], d, anchored['rigidity_anchor_state'])
    result = {'u_raw': d['u'].copy(), 'R_raw': d['R'].copy()}
    for name, model, x, y in [
        ('beta_only_frozen', unanchored['models']['beta_only_frozen'], xb, xr),
        ('joint_unanchored', unanchored['models']['joint_balanced'], xb, xr),
        ('joint_anchored', anchored['models']['anchored_joint_balanced'], ab, ar),
    ]:
        result['u_'+name] = d['u'] - x.dot(model['theta_beta'])
        result['R_'+name] = d['R'] * np.exp(-y.dot(model['theta_rigidity']))
    return result


def verify_replay(u, a):
    _, d = js.load_controlled_cache(BASE/'controlled_cache.npz')
    _, t = js.load_tracker_cache(HERE/'tracker_multifit_cache.npz')
    p = predict(d, t, u, a)
    with np.load(HERE/'joint_identifiability_predictions.npz') as old, np.load(HERE/'anchored_joint_predictions.npz') as anchored:
        for method, source, key in [('raw', old, 'raw'), ('beta_only_frozen', old, 'beta_only_frozen'),
                                    ('joint_unanchored', old, 'joint_balanced'),
                                    ('joint_anchored', anchored, 'anchored_joint_balanced')]:
            for field in ('u_', 'R_'):
                assert np.allclose(p[field+method], source[field+key], rtol=1e-12, atol=1e-12), (method, field)
    rc.log('Frozen feature/model replay matches all four original predictions.')


def evaluate(d, p):
    truth = d['truth_valid'] & np.isfinite(d['u_true']) & np.isfinite(d['R_true']) & (d['R_true']>0)
    rows = []
    for b in [-1, 0, 1, 2, 3]:
        for label in [-1, 0, 1]:
            mask = truth & ((d['bin']==b) if b>=0 else True) & ((d['label']==label) if label>=0 else True)
            for name in METHODS:
                u, r = p['u_'+name][mask], p['R_'+name][mask]
                labels = d['label'][mask]
                m2 = (6*r)**2*(u*u-1)
                vals = {'du':u-d['u_true'][mask], 'dbeta':1/u-1/d['u_true'][mask],
                        'dlogR':np.log(r/d['R_true'][mask]),
                        'mass':np.sqrt(m2[m2>=0]), 'mass2_error':m2-rc.MASSES[labels]**2}
                row = dict(method=name, bin=b, sample=['all','C12','C13'][label+1], n=int(mask.sum()),
                           mass2_auc=rc.auc(m2,labels) if label<0 else float('nan'),
                           rejection50=js.roc_rejection(m2,labels) if label<0 else float('nan'))
                for key, x in vals.items():
                    row['median_'+key] = float(np.median(x)); row['w68_'+key] = rc.w68(x)
                rows.append(row)
    rc.write_tsv(OUT/'metrics.tsv', rows)
    return rows


def main():
    started=time.time()
    OUT.mkdir(exist_ok=True)
    u = arrays(json.loads((HERE/'joint_weak_frozen_models.json').read_text()))
    a = arrays(json.loads((HERE/'anchored_joint_frozen_models.json').read_text()))
    verify_replay(u, a)
    used=set()
    for filename in ('manifest.tsv','blind_manifest.tsv'):
        with (BASE/filename).open() as f:
            used.update(row['path'] for row in csv.DictReader(f,delimiter='\t'))
    groups={}
    for label, pattern, count in [(0,rc.C12_PATTERN,256),(1,rc.C13_PATTERN,32)]:
        candidates=[p for p in sorted(glob.glob(pattern)) if p not in used]
        candidates.sort(key=lambda p: hashlib.sha256(('expanded-20260929:'+p).encode()).hexdigest())
        assert len(candidates)>=count
        groups[label]=candidates[:count]
    manifest=[dict(isotope=rc.ISOTOPE_NAMES[label],source=i,path=p) for label in (0,1) for i,p in enumerate(groups[label])]
    rc.write_tsv(OUT/'manifest.tsv',manifest)
    freeze={'model_sha256':{n:rc.sha256_file(HERE/n) for n in ['joint_weak_frozen_models.json','anchored_joint_frozen_models.json']},
            'script_sha256':rc.sha256_file(Path(__file__)), 'manifest_sha256':rc.sha256_file(OUT/'manifest.tsv'),
            'excluded_previous_files':len(used),'file_overlap':0,'selection':'256 C12 + 32 C13 unused files; deterministic SHA256 ordering; all passing events',
            'rigidity_edges_GV':rc.R_EDGES.tolist(),'training':'Original models and feature states frozen; no retraining or tuning.'}
    (OUT/'freeze.json').write_text(json.dumps(freeze,indent=2)+'\n')
    audit=[]
    pools=[rc.extract_pool(groups[label],label,[67,117][label],rc.MASSES[label],2,audit) for label in (0,1)]
    d=rc.concatenate(pools)
    del pools
    rc.write_tsv(OUT/'extraction_audit.tsv',audit)
    t, tracker_audit=rt.extract_tracker_for_split(d,{(2,label):groups[label] for label in (0,1)},'expanded')
    pieces=[]
    for first in range(0,len(d['u']),10000):
        part=slice(first,first+10000)
        pieces.append(predict(rc.subset(d,part),rc.subset(t,part),u,a))
    p=rc.concatenate(pieces)
    assert all(np.all(np.isfinite(v)) for v in p.values())
    d={k:d[k] for k in ['run','event','source','label','bin','R','u','truth_valid','u_true','R_true']}
    ids=np.rec.fromarrays([d[k] for k in ['label','run','event']])
    assert len(np.unique(ids))==len(ids), 'Duplicate expanded event IDs'
    for filename in ['predictions.npz','blind_predictions_with_truth.npz']:
        with np.load(BASE/filename) as old:
            oldids=np.rec.fromarrays([old[k] for k in ['label','run','event']])
            assert len(np.intersect1d(ids,oldids))==0, 'Event overlap with '+filename
    rc.write_tsv(OUT/'extraction_audit.tsv',audit)
    rc.write_tsv(OUT/'tracker_audit.tsv',tracker_audit)
    np.savez_compressed(OUT/'predictions.npz',**d,**p)
    rows=evaluate(d,p)
    result={'selected_events':len(d['u']),'truth_valid_events':int(d['truth_valid'].sum()),
            'source_files':len(manifest),'elapsed_seconds':time.time()-started,
            'predictions_sha256':rc.sha256_file(OUT/'predictions.npz'),
            'overall':[r for r in rows if r['bin']==-1 and r['sample']=='all']}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    rc.log('COMPLETED '+json.dumps(result))


if __name__=='__main__':
    main()
