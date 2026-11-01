#!/usr/bin/env python3
"""Measurement-only test of residual information after the scalar loss inverse."""
import json
from pathlib import Path
import numpy as np
import free_curve_timing as curve
import free_timing_transport as material
import joint_adc_likelihood as joint
import flight_conditional_likelihood as timing
import periodic_joint_adc as period
from research_amplitude_noise_20261013 import inputs,CONFIGS
from research_flight_conditional import write,digest
from research_conditional_mc import stamp

DEST=period.OUT/'curve_audit'


def stats(values):
    values=np.asarray(values,float);q=np.percentile(values,[0,2.5,16,50,84,97.5,100])
    return dict(mean=float(values.mean()),rms=float(np.sqrt(np.mean(values**2))),w68=float(.5*(q[4]-q[2])),percentiles=q)


def main():
    DEST.mkdir(exist_ok=True)
    if (DEST/'protocol.json').exists():raise FileExistsError('Curve audit was already started')
    write(DEST/'protocol.json',dict(created_utc=stamp(),source_sha256={name:digest(period.ROOT/'analysis'/name) for name in ['audit_free_curve_timing.py','free_curve_timing.py','free_timing_transport.py','joint_adc_likelihood.py']},
        hypothesis='The existing scalar material inverse sets d^T Q residual=0, whereas stationarity along the complete mean curve requires mu_prime(u)^T Q residual=0.',
        objective='Frozen seven-observable Mahalanobis residual using continuous measured-R mass and the same first-order mean material model.',
        estimator='Safeguarded Gauss-Newton refinement from the existing output; retain every event and all failure flags.',
        no_global_optimality_claim=True,no_beta_reference=True,no_model_retraining=True,
        scientific_scope='Numerical/output sensitivity audit, not a beta-resolution acceptance test. A nonlinear curve need not preserve exact finite noisy-signal translation response.'))
    rows=[]
    for dataset in ['flight','mc']:
        folder,a,g,adc,domain,train,test=inputs(dataset);select=np.flatnonzero(test)
        norm=json.loads((joint.OUT/dataset/'protocol.json').read_text())['feature_normalization']
        features,_=joint.response_features(a,train,norm);x=g['x'][select];r=a['R'][select];moments=g['moment'][select]
        for family,field,_ in CONFIGS[dataset]:
            source=joint.OUT/dataset/family/'freeze.json';theta=np.asarray(json.loads(source.read_text())['model']['theta'])
            y=joint.observation({field:a[field][select]},adc[select],field)
            rec=joint.reconstruct(theta,x,y,features[select]);q=joint.precision(theta)[1]
            weight=rec['weights'][:,:3].dot(timing.D)
            old=material.solve_center(rec['joint_u'],r,weight,moments,theta[:2])
            if not old['converged'].all():raise ValueError('Existing inverse did not converge')
            result=curve.refine(old['u'],rec['centered'],rec['direction'],q,r,moments,theta[:2])
            pick=np.arange(min(128,len(r)));v=result['u'][pick];delta=np.linspace(-.03,.03,len(pick));m0,m1=material.compress(moments[pick])
            displacement=curve.mean_curve(v+delta,r[pick],rec['direction'][pick],m0,m1,theta[:2])-curve.mean_curve(v,r[pick],rec['direction'][pick],m0,m1,theta[:2])
            changed=curve.refine(v+delta,rec['centered'][pick]+displacement,rec['direction'][pick],q,r[pick],moments[pick],theta[:2])
            response_error=changed['u']-v-delta
            path=DEST/(dataset+'_'+family+'.npz')
            np.savez_compressed(path,run=a['run'][select],event=a['event'][select],original_u=old['u'],**result)
            row=dict(dataset=dataset,family=family,events=len(r),converged=int(result['converged'].sum()),
                line_search_failures=int(result['line_search_failure'].sum()),iterations=result['iterations'],
                inverse_speed_change=stats(result['u']-old['u']),objective_reduction=stats(result['initial_objective']-result['objective']),
                maximum_stationarity_u=float(np.max(abs(result['score_in_u_units']))),
                noisy_mean_curve_injection_error=stats(response_error),prediction_sha256=digest(path),model_sha256=digest(source))
            rows.append(row);print('FREE_CURVE_AUDIT',json.dumps(row,default=lambda a:a.tolist()),flush=True)
    write(DEST/'evaluation.json',dict(complete=True,created_utc=stamp(),rows=rows,protocol_sha256=digest(DEST/'protocol.json'),
        no_beta_reference_read=True,production_replacement=False))


if __name__=='__main__':main()
