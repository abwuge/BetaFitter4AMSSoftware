#!/usr/bin/env python3
"""Refit quadrature only, preserving all statistical assumptions and old freezes."""
import argparse,json,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import joint_adc_likelihood as joint
import flight_conditional_likelihood as base
import free_timing_transport as transport
from joint_adc_scaled_optimizer import solve_scaled
from research_joint_adc import OUT as JOINT,ROOT,CONFIGS,inputs
from research_flight_conditional import load,write,digest
from research_conditional_mc import stamp

OUT=JOINT/'quadrature_refit'


def prepare():
    OUT.mkdir(exist_ok=True)
    if (OUT/'protocol.json').exists():raise FileExistsError('Quadrature refinement protocol exists')
    records=[]
    for dataset in ['mc','flight']:
        for name,field,nu in CONFIGS[dataset]:
            path=JOINT/dataset/name/'freeze.json';meta=json.loads((path.parent/'measurement_evaluation.json').read_text())
            if not meta['complete'] or digest(path)!=meta['freeze_sha256']:raise ValueError('Original joint model not frozen')
            records.append(dict(dataset=dataset,family=name,time_field=field,nu=nu,original_freeze=str(path),freeze_sha256=digest(path),protocol_sha256=digest(JOINT/dataset/'protocol.json')))
    write(OUT/'protocol.json',dict(created_utc=stamp(),records=records,nodes=[256,512],evaluation_nodes=[256,512,1024],
        rationale='Original 128-node event NLL integration errors reached a few hundredths; isolate quadrature-induced parameter/output changes',
        initialization='Sequential warm start from frozen 128-node model then complete 256-node fit; no beta truth or holdout metric choice',
        hypothesis_unchanged=['joint seven-dimensional Student response','same finite q domain','same ADC linear mean','same original training measurements','same bounds and penalty'],
        numerical_gate='Optimizer success and original-coordinate projected gradient max <1e-4; independently report dense-grid objective and output differences, no automatic production replacement',
        MC_truth_read=False,RICH_read=False,source_sha256=digest(Path(__file__)),optimizer_sha256=digest(ROOT/'analysis/joint_adc_scaled_optimizer.py')))
    print('JOINT_QUADRATURE_PROTOCOL',len(records),flush=True)


def one(record,protocol):
    dataset=record['dataset'];name=record['family'];field=record['time_field'];nu=record['nu']
    if digest(record['original_freeze'])!=record['freeze_sha256'] or digest(JOINT/dataset/'protocol.json')!=record['protocol_sha256']:raise ValueError('Original joint model changed')
    parent=JOINT/dataset/name;old=np.asarray(json.loads(Path(record['original_freeze']).read_text())['model']['theta']);p=json.loads((JOINT/dataset/'protocol.json').read_text())
    folder,a,g,adc,domain,train,test=inputs(dataset)
    if digest(folder/'measured.npz')!=p['measurement_sha256'] or digest(folder/'moments.npz')!=p['geometry_sha256']:raise ValueError('Quadrature measurements changed')
    features,norm=joint.response_features(a,train,p['feature_normalization']);theta=old;records=[];baseline=load(parent/'predictions_before_truth.npz')
    for nodes in protocol['nodes']:
        dest=OUT/dataset/name/str(nodes);dest.mkdir(parents=True,exist_ok=True)
        f=dest/'fit.json'
        if f.exists():result=json.loads(f.read_text())
        else:
            obj=joint.JointADC(a,g,train,features,adc,penalty=p['penalty'],nodes=nodes,time_field=field,rigidity_nu=nu,threads=4)
            result=solve_scaled(obj,theta,p['bounds']);write(f,result)
        print('JOINT_QUADRATURE_FIT',dataset,name,nodes,result['success'],result['penalized_nll'],result['projected_gradient_max'],result['elapsed_seconds'],flush=True)
        if not result['success'] or result['projected_gradient_max']>=1e-4:raise RuntimeError('Quadrature refit not converged')
        theta=np.asarray(result['theta']);ev=joint.JointADC(a,g,test,features,adc,penalty=0,nodes=nodes,time_field=field,rigidity_nu=nu,threads=4)
        pred=joint.reconstruct(theta,g['x'][test],ev.y,features[test]);out={k:a[k][test] for k in ['run','event','R']}
        for key in ['time_u','adc_u','joint_u','time_injection_gain']:out[key]=pred[key]
        for key,weight in [('time',pred['time_weights']),('joint',pred['weights'][:,:3].dot(base.D))]:
            inv=transport.solve_center(pred[key+'_u'],a['R'][test],weight,g['moment'][test],theta[:2])
            if not inv['converged'].all():raise ValueError('Quadrature loss inverse failed')
            out[key+'_transport_u']=inv['u']
        dense={};lastper=None
        for nq in protocol['evaluation_nodes']:
            evalobj=joint.JointADC(a,g,test,features,adc,penalty=0,nodes=nq,time_field=field,rigidity_nu=nu,threads=4)
            value,gradient,per=evalobj.evaluate(theta,events=True);out['nll_%d'%nq]=per
            dense[str(nq)]=dict(nll=value)
            if lastper is not None:dense[str(nq)]['difference_to_previous']=dict(rms=float(np.sqrt(np.mean((per-lastper)**2))),max_abs=float(np.max(abs(per-lastper))),mean=float(np.mean(per-lastper)))
            lastper=per
        np.savez_compressed(dest/'predictions.npz',**out)
        difference={}
        for key in ['joint_u','time_u','adc_u','time_transport_u','joint_transport_u']:
            delta=out[key]-baseline[key];difference[key]=dict(rms=float(np.sqrt(np.mean(delta*delta))),max_abs=float(np.max(abs(delta))),mean=float(delta.mean()),w68=float(np.diff(np.percentile(delta,[16,84]))[0]/2))
        row=dict(complete=True,created_utc=stamp(),nodes=nodes,events=int(test.sum()),fit_sha256=digest(f),output_sha256=digest(dest/'predictions.npz'),dense_evaluation=dense,
            output_difference_to_previous_nodes=difference,previous_nodes=128 if nodes==256 else 256,lambda_ADC=np.exp(theta[35:39]),material=theta[:2],MC_truth_read=False,RICH_read=False)
        write(dest/'evaluation.json',row);records.append(row);baseline=out
    final=dict(complete=True,created_utc=stamp(),dataset=dataset,family=name,records=records,protocol_sha256=digest(OUT/'protocol.json'))
    write(OUT/dataset/name/'complete.json',final);return final


def fit():
    p=json.loads((OUT/'protocol.json').read_text())
    if digest(Path(__file__))!=p['source_sha256'] or digest(ROOT/'analysis/joint_adc_scaled_optimizer.py')!=p['optimizer_sha256']:raise ValueError('Quadrature implementation changed')
    with ThreadPoolExecutor(max_workers=5) as pool:results=list(pool.map(lambda record:one(record,p),p['records']))
    write(OUT/'complete.json',dict(complete=True,created_utc=stamp(),records=results,protocol_sha256=digest(OUT/'protocol.json')))
    print('JOINT_QUADRATURE_COMPLETE',len(results),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','fit']);args=parser.parse_args()
    if args.stage=='prepare':prepare()
    else:fit()
