#!/usr/bin/env python3
"""Verify stored event identities, public inference, and transport conservation."""
import json
import subprocess
import sys
from pathlib import Path
import numpy as np
import rich_transport as rt
from rich_beta_estimator import RichBetaEstimator
import evaluate_rich_transport as ev


def main():
    dest=rt.ROOT/'results/rich_transport'
    tests=subprocess.run([sys.executable,'-m','unittest','discover','-s','tests','-p','test_rich_transport.py','-v'],
                         cwd=str(rt.ROOT),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    (dest/'unit_tests.txt').write_text(tests.stdout)
    if tests.returncode:raise RuntimeError(tests.stdout)
    result={'unit_test_returncode':tests.returncode,'datasets':{}}
    for kind in ['mc','data']:
        with np.load(str(dest/(kind+'_predictions.npz'))) as f:p={k:f[k] for k in f.files}
        with np.load(str(dest/(kind+'_reference_v2_masks.npz'))) as f:q={k:f[k] for k in f.files}
        assert np.array_equal(p['run'],q['run']) and np.array_equal(p['event'],q['event'])
        entry={'selected':len(p['run']),'RICH_v2':int(q['charge_consistent_reference'].sum())}
        for prefix in ['raw','new']:
            valid=p[prefix+'_valid'];loss=p[prefix+'_loss_GeV'][valid];planes=p[prefix+'_beta_planes'][valid]
            delta=p[prefix+'_initial_energy_GeV'][valid]-p[prefix+'_final_energy_GeV'][valid]
            err=float(np.max(np.abs(delta-loss.sum(axis=1))))
            increase=float(np.max(np.diff(planes,axis=1)))
            assert err<1e-10 and increase<1e-12 and np.all(loss>=0)
            entry[prefix+'_energy_conservation_max_abs_GeV']=err
            entry[prefix+'_max_beta_increase']=increase
        result['datasets'][kind]=entry
    with np.load('/tmp/betaEstimator-rich-transport/mc_evaluation_input.npz') as f:
        a={k:f[k][:16] for k in RichBetaEstimator.REQUIRED}
    model=RichBetaEstimator();first=model.predict(a,batch_size=16)
    a['rich_beta']=np.full((16,3),-.123);a['u_true']=np.full(16,123.);a['label']=np.arange(16)
    second=model.predict(a,batch_size=7)
    for key in first:np.testing.assert_allclose(first[key],second[key],rtol=0,atol=0,equal_nan=True)
    with np.load(str(dest/'mc_predictions.npz')) as f:
        dc=float(np.nanmax(np.abs(first['beta_center']-1/f['u_new'][:16])))
        dr=float(np.nanmax(np.abs(first['beta_rich_pred']-f['new_beta_planes'][:16,-1])))
        dl=float(np.nanmax(np.abs(first['loss_GeV']-f['new_loss_GeV'][:16])))
    assert max(dc,dr,dl)<1e-12
    result['public_API']={'events':16,'max_center_difference':dc,'max_RICH_difference':dr,'max_loss_difference_GeV':dl,
                          'ignored_RICH_truth_label_keys':True,'batch_size_invariance':True}
    audits={}
    for kind in ['mc','data']:
        rows=json.loads((dest/(kind+'_source_audit.json')).read_text())
        audits[kind]={'files':len(rows),'raw_entries':sum(r['raw'] for r in rows),'selected':sum(r['selected'] for r in rows)}
    result['source_counts']=audits
    result['code_sha256']={str(p.relative_to(rt.ROOT)):ev.digest(p) for p in list((rt.ROOT/'analysis').glob('*rich*.py'))+[rt.ROOT/'tests/test_rich_transport.py']}
    ev.dump(dest/'verification.json',result)
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
