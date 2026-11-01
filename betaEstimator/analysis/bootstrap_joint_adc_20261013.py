#!/usr/bin/env python3
"""Whole-run recalibration of prespecified primary joint-ADC configurations."""
import argparse,json,time,socket
from pathlib import Path
import numpy as np
import joint_adc_likelihood as joint
import flight_conditional_likelihood as base
import free_timing_transport as transport
from research_joint_adc import OUT as JOINT,ROOT,inputs
from research_additional_mc_20261013 import OUT as EXTRA
from research_flight_conditional import load,write,digest
from research_conditional_mc import stamp
from joint_adc_scaled_optimizer import solve_scaled

OUT=JOINT/'bootstrap'
CONFIGS=[('mc','native_studentR5','t_native',5.),('flight','newcalib_gaussianR','t',0.)]


def rebase(theta,old,new):
    theta=np.asarray(theta,float).copy();nf=(len(theta)-39)//4;alpha=theta[39:].reshape(4,nf)
    previous=alpha[:,1:6].copy();om=np.asarray(old['center']);os=np.asarray(old['scale']);nm=np.asarray(new['center']);ns=np.asarray(new['scale'])
    alpha[:,0]+=previous.dot((nm-om)/os);alpha[:,1:6]=previous*(ns/os)[None,:]
    return theta


def prepare():
    OUT.mkdir(exist_ok=True)
    if (OUT/'protocol.json').exists():raise FileExistsError('Joint bootstrap already declared')
    records=[]
    for dataset,name,field,nu in CONFIGS:
        qdir=JOINT/'quadrature_refit'/dataset/name
        if not json.loads((qdir/'complete.json').read_text())['complete']:raise ValueError('Joint quadrature check incomplete')
        check=json.loads((qdir/'512/evaluation.json').read_text())
        if check['output_difference_to_previous_nodes']['joint_transport_u']['rms']>1e-5 or check['output_difference_to_previous_nodes']['joint_transport_u']['max_abs']>1e-4:raise ValueError('Quadrature output not stable enough for bootstrap')
        folder,a,g,adc,domain,train,test=inputs(dataset);p=json.loads((JOINT/dataset/'protocol.json').read_text())
        records.append(dict(dataset=dataset,family=name,time_field=field,nu=nu,seed_path=str(qdir/'256/fit.json'),seed_sha256=digest(qdir/'256/fit.json'),
            original_protocol_sha256=digest(JOINT/dataset/'protocol.json'),train_runs=np.unique(a['run'][train]),train_events=int(train.sum()),test_events=int(test.sum()),
            measurement_sha256=digest(folder/'measured.npz'),geometry_sha256=digest(folder/'moments.npz'),quadrature_check_sha256=digest(qdir/'512/evaluation.json')))
    write(OUT/'protocol.json',dict(created_utc=stamp(),replicates=128,seed=926202610,nodes=256,records=records,
        selected_configs='MC Student rigidity, and independently fitted flight NewCalib Gaussian rigidity; primary configurations fixed, not selected by bootstrap or additional-run beta outcomes',
        sampling='Resample full measured-training run IDs with replacement; independent MC and flight draws; preserve every within-run event and duplicate selected runs',
        normalization='Recompute the five geometry normalization factors from each training resample; exact affine coefficient rebase before clipping, then complete refit',
        starts='Rebased dense-quadrature measured optimum and independent physical lambda3 plus resample-median ADC intercept; choose penalized training likelihood only',
        optimizer='Same curvature-scaled optimizer, both starts, success and original-coordinate projected gradient maximum<1e-4',
        outputs='Five free continuous speed branches and their extra-run predictions for MC; original evaluation events for flight. No isotope/source posterior in output.',
        uncertainty='These are genuine training-run recalibrations; later evaluation-run sampling is separate. Dense256 versus original128 numerical change is recorded, not hidden.',
        source_sha256=digest(Path(__file__)),optimizer_sha256=digest(ROOT/'analysis/joint_adc_scaled_optimizer.py'),
        extra_measurement_sha256=digest(EXTRA/'measured.npz'),extra_geometry_sha256=digest(EXTRA/'moments.npz'),
        extra_output_freeze_sha256=digest(EXTRA/'output_freeze.json'),MC_truth_read=False,RICH_read=False))
    print('JOINT_BOOTSTRAP_PROTOCOL',128,len(records),flush=True)


def prediction(theta,a,g,adc,features,select,field):
    small={field:a[field][select]};yy=joint.observation(small,adc[select],field);xx=g['x'][select];result=joint.reconstruct(theta,xx,yy,features[select])
    output={k:a[k][select] for k in ['run','event']}
    for k in ['joint_u','time_u','adc_u']:output[k]=result[k]
    for key,weight in [('time',result['time_weights']),('joint',result['weights'][:,:3].dot(base.D))]:
        inv=transport.solve_center(result[key+'_u'],a['R'][select],weight,g['moment'][select],theta[:2])
        if not inv['converged'].all():raise ValueError('Joint bootstrap material inverse failure')
        output[key+'_transport_u']=inv['u']
    if not all(np.isfinite(v).all() for k,v in output.items() if k.endswith('_u')):raise ValueError('Nonfinite bootstrap prediction')
    return output


def worker(rep):
    p=json.loads((OUT/'protocol.json').read_text())
    if not 0<=rep<p['replicates']:raise ValueError('Invalid joint bootstrap replicate')
    if digest(Path(__file__))!=p['source_sha256'] or digest(ROOT/'analysis/joint_adc_scaled_optimizer.py')!=p['optimizer_sha256']:raise ValueError('Joint bootstrap source changed')
    note=OUT/('%04d.json'%rep);output=note.with_suffix('.npz')
    if note.exists():
        old=json.loads(note.read_text())
        if old['complete'] and old['protocol_sha256']==digest(OUT/'protocol.json') and digest(output)==old['output_sha256']:return
        raise ValueError('Existing incomplete joint bootstrap')
    if digest(EXTRA/'measured.npz')!=p['extra_measurement_sha256'] or digest(EXTRA/'moments.npz')!=p['extra_geometry_sha256']:raise ValueError('Additional measured input changed')
    extra_a=load(EXTRA/'measured.npz');extra_g=load(EXTRA/'moments.npz');extra_select=np.ones(len(extra_a['run']),bool)
    start=time.monotonic();results=[];pred={}
    for ci,record in enumerate(p['records']):
        dataset=record['dataset'];name=record['family'];field=record['time_field'];nu=record['nu'];cp=json.loads((JOINT/dataset/'protocol.json').read_text())
        if digest(JOINT/dataset/'protocol.json')!=record['original_protocol_sha256'] or digest(record['seed_path'])!=record['seed_sha256']:raise ValueError('Joint bootstrap seed changed')
        folder,a,g,adc,domain,train,test=inputs(dataset)
        if digest(folder/'measured.npz')!=record['measurement_sha256'] or digest(folder/'moments.npz')!=record['geometry_sha256']:raise ValueError('Joint bootstrap data changed')
        rng=np.random.RandomState(p['seed']+rep+ci*100000);runs=np.asarray(record['train_runs']);draw=rng.choice(runs,len(runs),replace=True)
        index=np.concatenate([np.flatnonzero(train&(a['run']==r)) for r in draw]);features,norm=joint.response_features(a,index)
        old=np.asarray(json.loads(Path(record['seed_path']).read_text())['theta']);initial=rebase(old,cp['feature_normalization'],norm)
        bb=np.asarray(cp['bounds']);initial=np.clip(initial,bb[:,0],bb[:,1]);cold=joint.parameters_initial(features,adc,index,3.)
        obj=joint.JointADC(a,g,index,features,adc,penalty=cp['penalty'],nodes=p['nodes'],time_field=field,rigidity_nu=nu,threads=4)
        fits=[]
        for j,point in enumerate([initial,cold]):
            value=solve_scaled(obj,point,cp['bounds']);fits.append(value)
            print('JOINT_BOOTSTRAP_START',rep,dataset,name,j,value['success'],value['projected_gradient_max'],value['penalized_nll'],value['elapsed_seconds'],flush=True)
        best=min(fits,key=lambda row:row['penalized_nll']);theta=np.asarray(best['theta'])
        # Failed numerical fits are retained and marked, never hidden by beta outcome.
        converged=best['success'] and best['projected_gradient_max']<1e-4
        own=prediction(theta,a,g,adc,features,test,field)
        for k,v in own.items():pred[dataset+'_'+k]=v
        if dataset=='mc':
            eff,_=joint.response_features(extra_a,extra_select,norm);extra=prediction(theta,extra_a,extra_g,extra_a['adc'],eff,extra_select,'t_native')
            for k,v in extra.items():pred['extra_'+k]=v
        results.append(dict(dataset=dataset,family=name,drawn_runs=draw,training_events=len(index),normalization=norm,model=best,starts=fits,numerically_converged=converged))
    np.savez_compressed(output,**pred)
    write(note,dict(complete=True,created_utc=stamp(),replicate=rep,models=results,all_numerically_converged=all(r['numerically_converged'] for r in results),
        elapsed_seconds=time.monotonic()-start,host=socket.gethostname(),output_sha256=digest(output),protocol_sha256=digest(OUT/'protocol.json'),MC_truth_read=False,RICH_read=False))
    print('JOINT_BOOTSTRAP_COMPLETE',rep,time.monotonic()-start,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','worker']);parser.add_argument('--replicate',type=int);args=parser.parse_args()
    if args.stage=='prepare':prepare()
    else:worker(args.replicate)
