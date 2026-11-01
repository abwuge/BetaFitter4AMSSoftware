#!/usr/bin/env python3
"""Test quadrature domains and latent-state propagation without beta references."""
import argparse
import json
import time
from pathlib import Path
import numpy as np
import uproot
import flight_conditional_likelihood as model
import flight_mean_loss as physical
from research_flight_conditional import load,write,digest
from research_conditional_mc import stamp

ROOT=model.ROOT
BASE=ROOT/'results/flight_calibration_research/mean_loss_likelihood'
OUT=BASE/'domain_audit'
FAMILIES=[('native_gaussianR','t_native',0.),('newcalib_gaussianR','t',0.),('native_studentR5','t_native',5.)]


def stats(v):
    v=np.asarray(v,float)
    return dict(n=v.size,mean=float(v.mean()),rms=float(np.sqrt(np.mean(v*v))),
                percentiles=np.percentile(v,[0,2.5,16,50,84,97.5,100]).tolist())


def prepare():
    OUT.mkdir(exist_ok=True)
    if (OUT/'protocol.json').exists():raise FileExistsError('Domain audit protocol exists')
    a=load(BASE/'measured.npz');order=np.argsort(a['R']);ix=order[np.unique(np.linspace(0,len(order)-1,192).astype(int))]
    write(OUT/'protocol.json',dict(created_utc=stamp(),selected_rows=ix,selection='Uniform ranks in measured rigidity; no beta reference or time-residual selection',
        quadratures=[dict(name='original128',qrange=[.02,.7],nodes=128),dict(name='original512',qrange=[.02,.7],nodes=512),
                     dict(name='expanded512',qrange=[.005,1.4],nodes=512),dict(name='expanded1024',qrange=[.005,1.4],nodes=1024)],
        families=FAMILIES,latent_q_states=[.02,.1,.25,.4,.55,.7,1.,1.4],propagation_rows=ix[::12],
        material_scales=[[1.,1.],[2.,2.]],propagation_substeps=[1,4],
        model_freezes={name:digest(BASE/'conditional'/name/'freeze.json') for name,_,_ in FAMILIES},
        measured_sha256=digest(BASE/'measured.npz'),geometry_sha256=digest(BASE/'moments.npz'),
        scope='Frozen-parameter working-density sensitivity and numerical propagation, not a refit or detector-accuracy test',
        MC_truth_read=False,RICH_read=False))
    print('DOMAIN_PROTOCOL',len(ix),'events',flush=True)


def quadrature():
    protocol=json.loads((OUT/'protocol.json').read_text());a=load(BASE/'measured.npz');g=load(BASE/'moments.npz');ix=np.asarray(protocol['selected_rows'])
    if digest(BASE/'measured.npz')!=protocol['measured_sha256'] or digest(BASE/'moments.npz')!=protocol['geometry_sha256']:raise ValueError('Audit inputs changed')
    records=[];per={k:a[k][ix] for k in ['run','event','R']}
    for name,field,nu in FAMILIES:
        path=BASE/'conditional'/name/'freeze.json'
        if digest(path)!=protocol['model_freezes'][name]:raise ValueError('Frozen model changed')
        fr=json.loads(path.read_text());theta=np.asarray(fr['model']['theta']);models={}
        for rule in protocol['quadratures']:
            obj=model.ConditionalModel(a,g,ix,time_field=field,nodes=rule['nodes'],qrange=rule['qrange'],rigidity_nu=nu)
            loss,grad,nll=obj.evaluate(theta,events=True);models[rule['name']]=(obj,nll)
            per[name+'_'+rule['name']+'_nll']=nll
        original=models['original512'][1];expanded=models['expanded1024'][1];obj=models['expanded1024'][0]
        weights={k:[] for k in ['prior_below_original','prior_above_original','posterior_below_original','posterior_above_original','posterior_q_gt_04','posterior_q_gt_07','posterior_q_gt_10']}
        for i in range(len(ix)):
            result=obj.numpy_event(theta,i);pr=result['conditional_rigidity_weights'];po=result['posterior']
            for label,w in [('prior',pr),('posterior',po)]:
                weights[label+'_below_original'].append(w[:,obj.q<.02].sum());weights[label+'_above_original'].append(w[:,obj.q>.7].sum())
            for threshold,label in [(.4,'04'),(.7,'07'),(1.,'10')]:weights['posterior_q_gt_'+label].append(po[:,obj.q>threshold].sum())
        for key,value in weights.items():per[name+'_'+key]=np.array(value)
        record=dict(family=name,events=len(ix),original128_minus512=stats(models['original128'][1]-original),
            expanded512_minus1024=stats(models['expanded512'][1]-expanded),expanded_minus_original=stats(expanded-original),
            mass_q_weights={key:stats(v) for key,v in weights.items()})
        records.append(record);print('DOMAIN_QUADRATURE',name,'expanded-original',record['expanded_minus_original'],flush=True)
    np.savez_compressed(OUT/'quadrature_events.npz',**per)
    write(OUT/'quadrature.json',dict(complete=True,created_utc=stamp(),families=records,protocol_sha256=digest(OUT/'protocol.json'),events_sha256=digest(OUT/'quadrature_events.npz')))


def propagation():
    protocol=json.loads((OUT/'protocol.json').read_text());g=load(BASE/'moments.npz')
    with uproot.open(BASE/'geometry.root') as stream:raw=stream['geometry'].arrays(library='np')
    records=[];start=time.monotonic()
    for index in protocol['propagation_rows']:
        hi=np.asarray(raw['high'][index],float);lo=np.asarray(raw['low'][index],float);le=np.asarray(raw['length'][index],float)
        col=np.asarray(raw['abundance'][index],float).reshape(-1,9)*physical.ATOMIC_MASS
        group=np.isin(raw['segment'][index],[1,2]).astype(np.int64)
        args=(hi,lo,le,col,group,np.asarray(raw['at_high'][index],np.int64),np.asarray(raw['at_low'][index],np.int64),float(raw['center_z'][index]),raw['tof_z'][index].astype(float))
        x=g['x'][index];xc=x-x.mean();w=xc/np.dot(xc,xc)
        for mass in model.MASSES:
            for q in protocol['latent_q_states']:
                u=float(np.hypot(1,mass*q/6))
                for kappa in protocol['material_scales']:
                    fine,loss=physical.propagate_cells(u,mass,kappa,*args,substeps=4)
                    coarse,_=physical.propagate_cells(u,mass,kappa,*args,substeps=1)
                    linear=u*x+physical.first_order(u,mass,g['moment'][index]).dot(kappa)
                    valid=bool(np.isfinite(fine).all());row=dict(row=int(index),mass=float(mass),q=q,kappa=kappa,valid=valid,
                        u=u,linear_slope=float(w.dot(linear)),linear_max_contrast=float(np.max(abs(model.D.dot(linear-u*x)))))
                    if valid:row.update(linear_minus_nonlinear_u=float(w.dot(linear-fine)),max_time_contrast_error_ns=float(np.max(abs(model.D.dot(linear-fine)))),
                        one_minus_four_substep_u=float(w.dot(coarse-fine)),loss_GeV=loss)
                    records.append(row)
        print('DOMAIN_PROPAGATION',index,'configurations',len(records),'seconds',time.monotonic()-start,flush=True)
    summary=[]
    for q in protocol['latent_q_states']:
        for scales in protocol['material_scales']:
            rows=[r for r in records if r['q']==q and r['kappa']==scales];valid=[r for r in rows if r['valid']]
            summary.append(dict(q=q,kappa=scales,configurations=len(rows),valid=len(valid),stopped=len(rows)-len(valid),
                linear_minus_nonlinear_u=stats([r['linear_minus_nonlinear_u'] for r in valid]) if valid else None,
                time_contrast_error_ns=stats([r['max_time_contrast_error_ns'] for r in valid]) if valid else None,
                one_minus_four_substep_u=stats([r['one_minus_four_substep_u'] for r in valid]) if valid else None))
    write(OUT/'propagation.json',dict(complete=True,created_utc=stamp(),summary=summary,records=records,elapsed_seconds=time.monotonic()-start,
        protocol_sha256=digest(OUT/'protocol.json'),scope='Same elemental Bethe and material map; tests numerical linearization only'))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','quadrature','propagation']);args=parser.parse_args()
    if args.stage=='prepare':prepare()
    elif args.stage=='quadrature':quadrature()
    else:propagation()
