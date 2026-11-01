#!/usr/bin/env python3
"""Whole-training-run refits of the frozen MC model families."""
import argparse
import json
import time
from pathlib import Path
import numpy as np
import flight_conditional_likelihood as f
import free_timing_transport as transport
from research_conditional_mc import OUT,ROOT,inputs,stamp
from research_conditional_mc_reference import check_freeze
from research_flight_conditional import write,digest,optimize

BOOT=OUT/'bootstrap'
LABELS=[('native_gaussianR',0),('native_studentR5',5)]
VARIANTS=[('fitted','model',{}),('zero_loss','null_model',{0:0.,1:0.}),('unit_material','unit_material_model',{0:1.,1:1.})]


def prepare(count):
    check_freeze();BOOT.mkdir(exist_ok=True)
    if (BOOT/'protocol.json').exists():raise FileExistsError('Bootstrap protocol already frozen')
    a,g=inputs();train=a['fold']!=0;test=~train
    value=dict(created_utc=stamp(),replicates=count,seed=720261011,sampling='Training run IDs sampled with replacement; all events of each drawn run retained, including duplicate runs',
        family_labels=[x[0] for x in LABELS],variants=[x[0] for x in VARIANTS],
        initialization='Free-material fit uses both the original MC optimum and an independent zero-material 50/50-mixture start; nuisance controls use their own original MC optimum',
        selection='Minimum resampled training NLL only; no evaluation truth or reference files opened by workers',
        train_runs=np.unique(a['run'][train]),test_runs=np.unique(a['run'][test]),nodes=128,
        measured_sha256=digest(OUT/'measured.npz'),geometry_sha256=digest(OUT/'moments.npz'),
        output_freeze_sha256=digest(OUT/'output_freeze.json'),script_sha256=digest(Path(__file__)),
        uncertainties='Training-run refit variability; evaluation-run resampling will be added independently at final scoring; detector structure uncertainty is separate',
        RICH_read=False,MC_truth_read=False)
    write(BOOT/'protocol.json',value);print('BOOTSTRAP_PROTOCOL',count,flush=True)


def replicate(rep,a,g,protocol):
    note=BOOT/('%04d.json'%rep);output=BOOT/('%04d.npz'%rep)
    if note.exists():
        old=json.loads(note.read_text())
        if old.get('complete') and digest(output)==old['output_sha256']:return
        raise ValueError('Incomplete bootstrap result already exists')
    runs=np.asarray(protocol['train_runs']);rng=np.random.RandomState(protocol['seed']+rep)
    draw=rng.choice(runs,len(runs),replace=True)
    idx=np.concatenate([np.flatnonzero((a['run']==run)&(a['fold']!=0)) for run in draw])
    test=a['fold']==0;x=g['x'][test];t=a['t_native'][test];R=a['R'][test];moment=g['moment'][test]
    results=[];pred={k:a[k][test] for k in ['run','event']}
    start=time.monotonic()
    for label,nu in LABELS:
        frozen=json.loads((OUT/label/'freeze.json').read_text())
        model=f.ConditionalModel(a,g,idx,nodes=protocol['nodes'],rigidity_nu=nu,threads=4)
        for variant,key,fixed in VARIANTS:
            initial=np.asarray(frozen[key]['theta']);fits=[optimize(model,initial,fixed=fixed,label='boot_%d_%s_%s'%(rep,label,variant))]
            if not fixed:fits.append(optimize(model,f.initial(0,.5),label='boot_%d_%s_independent'%(rep,label)))
            best=min(fits,key=lambda v:v['nll'])
            if not best['success'] or best['projected_gradient_max']>1e-3:
                retry=optimize(model,best['theta'],fixed=fixed,label='boot_retry',maxiter=1200)
                fits.append(retry);best=min(fits,key=lambda v:v['nll'])
            theta=np.asarray(best['theta']);w=transport.weights(theta,x);s=np.sum(w*t,axis=1)
            center=transport.solve_center(s,R,w,moment,theta[:2])
            if not center['converged'].all():raise ValueError('Bootstrap transport inverse failed')
            tag=label+'_'+variant;pred[tag+'_u']=s;pred[tag+'_transport_u']=center['u']
            results.append(dict(label=tag,best=best,starts=fits,max_transport_inverse_residual=float(np.max(np.abs(center['residual'])))))
    np.savez_compressed(output,**pred)
    write(note,dict(complete=True,replicate=rep,created_utc=stamp(),drawn_runs=draw,training_events=len(idx),
        models=results,elapsed_seconds=time.monotonic()-start,output_sha256=digest(output),protocol_sha256=digest(BOOT/'protocol.json'),
        RICH_read=False,MC_truth_read=False))
    print('BOOTSTRAP_DONE',rep,'seconds',time.monotonic()-start,flush=True)


def worker(index,workers):
    check_freeze();protocol=json.loads((BOOT/'protocol.json').read_text())
    if digest(Path(__file__))!=protocol['script_sha256']:raise ValueError('Bootstrap source changed')
    if digest(OUT/'measured.npz')!=protocol['measured_sha256'] or digest(OUT/'moments.npz')!=protocol['geometry_sha256']:raise ValueError('Bootstrap measured inputs changed')
    a,g=inputs()
    for rep in range(index,protocol['replicates'],workers):replicate(rep,a,g,protocol)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','worker']);p.add_argument('--replicates',type=int,default=128)
    p.add_argument('--worker',type=int,default=0);p.add_argument('--workers',type=int,default=50);args=p.parse_args()
    if args.stage=='prepare':prepare(args.replicates)
    else:worker(args.worker,args.workers)
