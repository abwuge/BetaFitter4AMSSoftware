#!/usr/bin/env python3
"""Refit the conditional time model on independently stored AMS MC measurements.

Measured preparation and fitting are separate from reference extraction. The
existing flight models are read only for opaque integrity hashes.
"""
import argparse
import csv
import hashlib
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import uproot
import flight_conditional_likelihood as f
import research_flight_material as material
from research_flight_conditional import load, write, digest, optimize
from research_conditional_output import weights

ROOT=f.ROOT
OUT=ROOT/'results/flight_calibration_research/conditional_mc'
MANIFEST=ROOT/'results/c_joint_beta_rigidity_closure/expanded_validation/manifest.tsv'
FIELDS=['run','event','itrtrack','ibetah','tk_rigidityk','tk_erigidity1',
        'tk_rigidity1','tk_chis1','tk_q','tk_qin','tk_pos','tk_dir',
        'tof_tl','tof_leng','tof_pos','tof_barid','tof_hsumh','tof_hsumhu',
        'tof_bpatt','tof_edep']
OPTIONAL=['tof_etl','tof_adca','tof_betah']
C=29.9792458


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())


def measured_record(record):
    with uproot.open(record['path']) as stream:
        tree=stream['amstreea'];available=set(tree.keys())
        missing=set(FIELDS)-available
        if missing:raise ValueError('Missing measured fields: %s %s'%(record['path'],missing))
        fields=FIELDS+[k for k in OPTIONAL if k in available]
        a=tree.arrays(fields,entry_stop=record['entry_stop'],library='np')
    n=len(a['run']);R=a['tk_rigidityk'][:,1,1].astype(float)
    sigma=np.abs(a['tk_erigidity1'][:,2,1].astype(float));rel=sigma*np.abs(R)
    charge=np.where(a['tk_qin'][:,0,2]<2.5,a['tk_q'][:,1],a['tk_qin'][:,0,2])
    t=a['tof_tl'].astype(float);x=a['tof_leng'].astype(float)/C
    xc=x-x.mean(1)[:,None]
    with np.errstate(invalid='ignore',divide='ignore'):
        u=np.sum(xc*t,axis=1)/np.sum(xc*xc,axis=1)
    chi=a['tk_chis1'][:,2,2,1,:]
    good=np.isfinite(charge)&(charge>=5.5)&(charge<6.5)
    good&=np.isfinite(R)&(R>=2.5)&(R<6)&np.isfinite(rel)&(rel>0)&(rel<.2)
    good&=np.isfinite(chi).all(1)&((chi>0)&(chi<=10)).all(1)
    good&=(a['tof_hsumh']==4)&(a['tof_hsumhu']==4)&(a['tof_bpatt']==4444)
    good&=np.isfinite(t).all(1)&(np.abs(t)<1000).all(1)&np.isfinite(x).all(1)&(np.ptp(x,axis=1)>1)
    good&=np.isfinite(u)&(u>0)&((a['tof_barid']>=0)&(a['tof_barid']<10)).all(1)
    good&=np.isfinite(a['tk_pos']).all((1,2))&np.isfinite(a['tk_dir']).all((1,2))
    good&=(a['tk_pos'][:,4,2]>0)&(a['tk_pos'][:,5,2]<0)&(np.abs(a['tk_pos'][:,4,:])<200).all(1)
    good&=(np.abs(a['tk_dir'][:,4,2])>.1)&(a['itrtrack']>=0)&(a['ibetah']>=0)
    good&=np.isfinite(a['tof_edep']).all(1)&(a['tof_edep']>0).all(1)
    data=dict(run=a['run'],event=a['event'],entry=np.arange(n),source=np.full(n,record['source'],int),
              itrtrack=a['itrtrack'],ibetah=a['ibetah'],R=R,relR=rel,
              R_default=a['tk_rigidity1'][:,2,2,1],charge=charge,
              t_native=t,t=t.copy(),x=x,bar=a['tof_barid'],tk_pos=a['tk_pos'],tk_dir=a['tk_dir'],
              tof_pos=a['tof_pos'],edep=a['tof_edep'],chi=chi)
    for field,shape in [('tof_etl',(4,)),('tof_adca',(4,2)),('tof_betah',())]:
        data[field]=a.get(field,np.full((n,)+shape,np.nan))
    data={k:v[good] for k,v in data.items()}
    note=dict(source=record['source'],input_events=n,selected_events=int(good.sum()),
              fields_read=fields,optional_missing=[k for k in OPTIONAL if k not in fields])
    return data,note


def prepare(entry_stop,limit,workers):
    OUT.mkdir(exist_ok=True)
    if (OUT/'measured.npz').exists():raise FileExistsError('MC measured snapshot already exists')
    rows=list(csv.DictReader(MANIFEST.open(),delimiter='\t'))
    records=[dict(source=i,path=r['path'],entry_stop=entry_stop) for i,r in enumerate(rows)]
    protocol=dict(created_utc=stamp(),manifest_sha256=digest(MANIFEST),records=records,
        measured_fields=FIELDS,optional_fields=OPTIONAL,rigidity_GV=[2.5,6],maximum_events=limit,
        selection='Measured carbon charge, Kalman V6 inner relative error <0.2 and chi-square <=10, native four fired/used good-pattern times, valid geometry',
        sample='Pre-existing C12/C13 MC development source manifest; no per-event generated isotope or truth used for preparation or fitting',
        split='SHA256 of decimal run modulo5; fold0 evaluation, folds1-4 training',
        rigidity='tk_rigidityk[1,1] = Kalman V6 inner GetRigidity(fit,1), native z=0 convention',
        rigidity_error='abs(tk_erigidity1[2,1]) = same Kalman V6 inner fit GetErrRinv; reference-plane response requires audit',
        geometry='Stored native Tracker states and TOF path from the same event association; BuildProfile reference remains L5 z=1.698 cm',
        reference_plan='Audit native primary MCEventg acceptance records and their GTrackID, parentID, particle and mass before L5 speed evaluation',
        reference_plane_caveat='Measured R is at z=0 while material time model is referenced to L5; keep this approximation explicit',
        fit_plan=['native_gaussianR','native_studentR5'],initializations=[[1,.1],[0,.5],[2,.9]],
        fixed_controls=[[0,0],[1,1]],output='Free GLS time slope; no event isotope posterior; constant u=1 negative control',
        truth_read=False,RICH_read=False,flight_coefficients_imported=False)
    write(OUT/'protocol.json',protocol)
    pieces=[];audit=[]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for i,(a,note) in enumerate(pool.map(measured_record,records),1):
            pieces.append(a);audit.append(note)
            if i%24==0 or i==len(records):print('MEASURED_MC',i,'/',len(records),'selected',sum(x['selected_events'] for x in audit),flush=True)
    a={k:np.concatenate([p[k] for p in pieces]) for k in pieces[0]}
    ids=material.keys(a)
    if len(np.unique(ids))!=len(ids):raise ValueError('Duplicate run/event in MC sources')
    parent=len(ids)
    ranks=np.array([int.from_bytes(hashlib.sha256(('%d:%d'%(rr,ee)).encode()).digest()[:8],'big') for rr,ee in zip(a['run'],a['event'])],np.uint64)
    ix=np.sort(np.argsort(ranks)[:min(limit,parent)]);a={k:v[ix] for k,v in a.items()}
    lookup={int(rr):int(hashlib.sha256(str(int(rr)).encode()).hexdigest(),16)%5 for rr in np.unique(a['run'])}
    a['fold']=np.array([lookup[int(rr)] for rr in a['run']],np.int8)
    np.savez_compressed(OUT/'measured.npz',**a)
    values=np.column_stack([a['run'],a['event'],a['R'],a['bar'],a['x']*C,a['tk_pos'].reshape(len(ix),-1),a['tk_dir'].reshape(len(ix),-1)])
    np.savetxt(OUT/'geometry_input.txt',values,fmt=['%d','%d','%.12g']+['%d']*4+['%.12g']*58)
    write(OUT/'preparation.json',dict(complete=True,created_utc=stamp(),events=len(ix),parent_selected=parent,
        train_events=int((a['fold']!=0).sum()),test_events=int((a['fold']==0).sum()),
        train_runs=len(np.unique(a['run'][a['fold']!=0])),test_runs=len(np.unique(a['run'][a['fold']==0])),
        records=audit,measured_sha256=digest(OUT/'measured.npz'),geometry_input_sha256=digest(OUT/'geometry_input.txt'),
        isotope_labels_stored_in_measured=False,truth_read=False,RICH_read=False))
    print('PREPARED_MC',len(ix),flush=True)


def inputs():
    a=load(OUT/'measured.npz');g=load(OUT/'moments.npz')
    if not np.array_equal(material.keys(a),material.keys(g)):raise ValueError('MC geometry identity mismatch')
    if not g['valid'].all():raise ValueError('Invalid MC geometry')
    return a,g


def fit_model(nu,nodes,starts):
    a,g=inputs();train=a['fold']!=0;test=~train
    if np.intersect1d(a['run'][train],a['run'][test]).size:raise ValueError('Run overlap')
    label='native_gaussianR' if nu<=0 else 'native_studentR5';folder=OUT/label;folder.mkdir(exist_ok=True)
    if (folder/'freeze.json').exists():raise FileExistsError('MC fit already frozen')
    model=f.ConditionalModel(a,g,train,nodes=nodes,rigidity_nu=nu)
    write(folder/'protocol.json',dict(created_utc=stamp(),parameters=f.PARAMETERS,bounds=f.BOUNDS,nodes=nodes,
        qrange=model.qrange,rigidity_nu=nu,time_nu=5,starts=starts,train_runs=np.unique(a['run'][train]),test_runs=np.unique(a['run'][test]),
        measured_sha256=digest(OUT/'measured.npz'),geometry_sha256=digest(OUT/'moments.npz'),
        model_module_sha256=digest(Path(f.__file__)),library_sha256=digest(f.OUT/'libconditional.so'),
        truth_read=False,RICH_read=False,flight_coefficients_imported=False))
    results=[]
    for j,(k,p) in enumerate([(1,.1),(0,.5),(2,.9)][:starts]):
        result=optimize(model,f.initial(k,p),label=label+'_mc_start%d'%j)
        results.append(result);write(folder/('start_%d.json'%j),result)
        print('MC_START_DONE',label,j,result['nll'],result['theta'][:2],result['projected_gradient_max'],flush=True)
    best=min(results,key=lambda z:z['nll'])
    null=optimize(model,best['theta'],fixed={0:0.,1:0.},label=label+'_mc_null')
    unit=optimize(model,best['theta'],fixed={0:1.,1:1.},label=label+'_mc_unit')
    freeze=dict(created_utc=stamp(),model=best,null_model=null,unit_material_model=unit,starts=results,
        selection='Minimum training NLL across prescribed starts',protocol_sha256=digest(folder/'protocol.json'),truth_read=False,
        RICH_read=False,flight_coefficients_imported=False)
    write(folder/'freeze.json',freeze)
    pred={k:a[k][test] for k in ['run','event','source','entry']};evaluation={}
    x=g['x'][test];t=a['t_native'][test];xc=x-x.mean(1)[:,None];w0=xc/np.sum(xc*xc,axis=1)[:,None]
    pred['ols_u']=np.sum(w0*t,axis=1);pred['constant_u']=np.ones(len(x))
    pred['native_beta']=a['tof_betah'][test];inject={}
    for name,key in [('fitted','model'),('zero_loss','null_model'),('unit_material','unit_material_model')]:
        theta=np.array(freeze[key]['theta']);w=weights(theta,x);u=np.sum(w*t,axis=1);pred[name+'_u']=u
        aa=np.linspace(-2,2,len(t));bb=np.linspace(-.2,.2,len(t))
        err=np.sum(w*(t+aa[:,None]+bb[:,None]*x),axis=1)-u-bb
        inject[name]=dict(max_error=float(np.max(np.abs(err))),median_weights=np.median(w,axis=0))
        if np.max(np.abs(err))>1e-10:raise ValueError('MC coherent injection failed')
    for nq in [nodes,2*nodes]:
        m=f.ConditionalModel(a,g,test,nodes=nq,rigidity_nu=nu);evaluation[str(nq)]={}
        for name,key in [('fitted','model'),('zero_loss','null_model'),('unit_material','unit_material_model')]:
            loss,grad,per=m.evaluate(np.asarray(freeze[key]['theta']),events=True)
            evaluation[str(nq)][name]=loss;pred['nll_%s_%d'%(name,nq)]=per
    np.savez_compressed(folder/'predictions_before_truth.npz',**pred)
    write(folder/'measured_evaluation.json',dict(complete=True,created_utc=stamp(),events=int(test.sum()),runs=len(np.unique(a['run'][test])),
        mean_conditional_nll=evaluation,injection=inject,freeze_sha256=digest(folder/'freeze.json'),
        predictions_sha256=digest(folder/'predictions_before_truth.npz'),truth_read=False))
    print('MC_FROZEN',label,evaluation,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','moments','fit']);p.add_argument('--entry-stop',type=int,default=8000)
    p.add_argument('--limit',type=int,default=24000);p.add_argument('--workers',type=int,default=4)
    p.add_argument('--nu',type=float,default=0);p.add_argument('--nodes',type=int,default=128);p.add_argument('--starts',type=int,default=3)
    args=p.parse_args()
    if args.stage=='prepare':prepare(args.entry_stop,args.limit,args.workers)
    elif args.stage=='moments':material.OUT=OUT;material.moments()
    else:fit_model(args.nu,args.nodes,args.starts)
