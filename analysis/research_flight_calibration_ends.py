#!/usr/bin/env python3
"""Flight-only, out-of-run test of two-ended relative amplitude response.

Uses no beta, charge, RICH, or MC field. This diagnoses a response contrast,
not absolute deposited energy or beta resolution. Position is extrapolated
from the nearest stored Tracker state, not inferred from TOF times.
"""
import json
import hashlib
from pathlib import Path
import numpy as np
import uproot

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/flight_calibration_research'
FIELDS=['run','event','tof_adca','tof_barid','tof_pos','tk_pos','tk_dir']
PROTOCOL=dict(fields=FIELDS,source_count=9,entries_per_file=5000,
    split='SHA256(str(run)) modulo 3 equals 0 for evaluation; otherwise fit',
    shared_model='Per bar: intercept plus quadratic in a training-only principal coordinate',
    time_model='Same position coefficients plus source-epoch intercepts; first epoch reference',
    minimum_train_per_epoch=30,minimum_evaluation_per_epoch=10,
    robust_iterations=20,truth_used=False,RICH_used=False,beta_used=False,
    caveats='All observed charges; no established saturation or trigger-efficiency correction. First-5000-entry sampling is not a representative flux sample. ADC ratio cancellation is approximate outside linear equal-response operation.')


def robust_fit(x,y):
    beta=np.linalg.lstsq(x,y,rcond=1e-10)[0]
    for _ in range(PROTOCOL['robust_iterations']):
        residual=y-x.dot(beta)
        scale=max(1e-5,1.4826*np.median(np.abs(residual-np.median(residual))))
        weight=np.minimum(1.,2.5*scale/np.maximum(np.abs(residual),1e-12))
        revised=np.linalg.lstsq(x*np.sqrt(weight)[:,None],y*np.sqrt(weight),rcond=1e-10)[0]
        if np.max(np.abs(revised-beta))<1e-9:beta=revised;break
        beta=revised
    return beta


def stats(v):
    q=np.percentile(v,[16,50,84])
    return dict(n=len(v),mean=float(v.mean()),median=float(q[1]),w68=float((q[2]-q[0])/2),rmse=float(np.sqrt(np.mean(v*v))))


def read_source(record,epoch):
    with uproot.open(record['path']) as f:
        tree=f['amstreea'];a=tree.arrays(FIELDS,entry_stop=min(5000,tree.num_entries),library='np')
    n=len(a['run']);pieces=[]
    pos=a['tk_pos'].astype(float);direction=a['tk_dir'].astype(float)
    tracker_ok=np.isfinite(pos).all(axis=2)&np.isfinite(direction).all(axis=2)&(np.abs(direction[:,:,2])>.05)
    tracker_ok&=(np.max(np.abs(pos),axis=2)<2000)&(np.linalg.norm(direction,axis=2)>.1)
    for layer in range(4):
        z=a['tof_pos'][:,layer,2].astype(float)
        distance=np.where(tracker_ok,np.abs(pos[:,:,2]-z[:,None]),np.inf)
        nearest=np.argmin(distance,axis=1)
        p=pos[np.arange(n),nearest];d=direction[np.arange(n),nearest]
        with np.errstate(invalid='ignore',divide='ignore'):
            xy=p[:,:2]+(z-p[:,2])[:,None]*d[:,:2]/d[:,2,None]
            adc=a['tof_adca'][:,layer,:].astype(float)
            y=np.log(adc[:,0]/adc[:,1])
        bar=a['tof_barid'][:,layer]
        ok=(bar>=0)&(bar<10)&np.isfinite(adc).all(axis=1)&(adc>0).all(axis=1)
        ok&=np.isfinite(xy).all(axis=1)&(np.abs(xy)<150).all(axis=1)&np.isfinite(y)
        ok&=(np.min(distance,axis=1)<80)&np.isfinite(z)&(np.abs(z)<100)
        pieces.append(dict(run=a['run'][ok],event=a['event'][ok],bar=bar[ok].astype(int),
            layer=np.full(ok.sum(),layer,int),epoch=np.full(ok.sum(),epoch,int),xy=xy[ok],y=y[ok]))
    return {k:np.concatenate([p[k] for p in pieces]) for k in pieces[0]}


def main():
    inventory=json.loads((OUT/'input_inventory.json').read_text())
    (OUT/'ends_protocol.json').write_text(json.dumps(PROTOCOL,indent=2)+'\n')
    pieces=[read_source(record,i) for i,record in enumerate(inventory['records'])]
    data={k:np.concatenate([p[k] for p in pieces]) for k in pieces[0]}
    ids=np.rec.fromarrays([data['run'],data['event'],data['layer']],names='run,event,layer')
    if len(np.unique(ids))!=len(ids):raise ValueError('Duplicate event-layer identity')
    allruns=np.unique(data['run'])
    holdout={int(r) for r in allruns if int(hashlib.sha256(str(int(r)).encode()).hexdigest(),16)%3==0}
    evaluation=np.isin(data['run'],list(holdout))
    assert not set(data['run'][~evaluation]).intersection(data['run'][evaluation])
    models=[];pred0=np.full(len(evaluation),np.nan);pred1=pred0.copy()
    for layer in range(4):
        for bar in range(10):
            channel=(data['layer']==layer)&(data['bar']==bar)
            supported=[]
            for epoch in range(9):
                here=channel&(data['epoch']==epoch)
                if (here&~evaluation).sum()>=30 and (here&evaluation).sum()>=10:supported.append(epoch)
            if len(supported)<3:continue
            select=channel&np.isin(data['epoch'],supported)
            train=select&~evaluation;test=select&evaluation
            center=data['xy'][train].mean(axis=0)
            _,_,vt=np.linalg.svd(data['xy'][train]-center,full_matrices=False)
            axis=vt[0]
            coordinate=(data['xy']-center).dot(axis)/100
            base=np.column_stack([np.ones(len(coordinate)),coordinate,coordinate*coordinate])
            design=np.column_stack([base]+[(data['epoch']==e).astype(float) for e in supported[1:]])
            if np.linalg.matrix_rank(design[train])!=design.shape[1]:continue
            b0=robust_fit(base[train],data['y'][train]);b1=robust_fit(design[train],data['y'][train])
            pred0[test]=base[test].dot(b0);pred1[test]=design[test].dot(b1)
            models.append(dict(layer=layer,bar=bar,epochs=supported,train_n=int(train.sum()),evaluation_n=int(test.sum()),
                center=center.tolist(),coordinate_axis=axis.tolist(),pooled=b0.tolist(),time_conditioned=b1.tolist()))
    # Save all coefficients before calculating held-out residual summaries.
    freeze=dict(protocol=PROTOCOL,models=models,training_runs=sorted(set(int(x) for x in data['run'][~evaluation])),
        evaluation_runs=sorted(holdout),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'ends_freeze.json').write_text(json.dumps(freeze,indent=2)+'\n')
    valid=evaluation&np.isfinite(pred0)&np.isfinite(pred1)
    rows=[]
    for layer in range(4):
        m=valid&(data['layer']==layer)
        if not m.sum():continue
        rows.append(dict(layer=layer,pooled=stats(data['y'][m]-pred0[m]),time_conditioned=stats(data['y'][m]-pred1[m])))
    result=dict(complete=True,measured_hit_rows=len(evaluation),eligible_evaluation_hit_rows=int(evaluation.sum()),
        evaluated_common_hit_rows=int(valid.sum()),model_count=len(models),rows=rows,
        training_run_count=len(freeze['training_runs']),evaluation_run_count=len(holdout),run_overlap=0,
        RICH_used=False,MC_used=False,beta_used=False,absolute_energy_calibration=False,beta_resolution_evaluation=False,
        metric='Log ADC end-ratio residual on identical out-of-run hits',
        limits=PROTOCOL['caveats']+' Nearest-state straight extrapolation and quadratic position response are diagnostic approximations; epoch differences are not uniquely attributed to gain aging.',
        script_sha256=freeze['script_sha256'],freeze_sha256=hashlib.sha256((OUT/'ends_freeze.json').read_bytes()).hexdigest())
    (OUT/'ends_response.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
