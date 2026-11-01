#!/usr/bin/env python3
"""Two-scale projected-time diagnostic, not the full physical likelihood.

Measured rigidity defines design states under two separately reported mass
hypotheses. It is not treated as a truth label, and its migration is NOT
modeled here. Parameters from this diagnostic are not production calibration.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
import flight_mean_loss as m
from flight_timing_redundancy import projector

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/flight_calibration_research/mean_loss_likelihood'

def native(v):
    if isinstance(v,dict):return {str(k):native(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [native(x) for x in v]
    if isinstance(v,np.ndarray):return native(v.tolist())
    if isinstance(v,np.generic):return v.item()
    return v

def stats(a):
    a=np.asarray(a,float);q=np.percentile(a,[16,50,84]);return dict(n=len(a),mean=float(a.mean()),median=float(q[1]),w68=float((q[2]-q[0])/2),rmse=float(np.sqrt(np.mean(a*a))))

def fit(b,y,scale,base_weight=None,bounds=None):
    n=len(y);base=np.ones(n) if base_weight is None else np.asarray(base_weight,float)
    weight=base.copy();coef=np.zeros(2);h=np.einsum('nli,nlj->nij',b,b);v=np.einsum('nli,nl->ni',b,y)
    for iteration in range(100):
        f=np.einsum('n,nij->ij',weight,h);rhs=np.einsum('n,ni->i',weight,v)
        candidate=np.linalg.solve(f,rhs)
        if bounds is not None:
            lo,hi=bounds;solutions=[]
            if np.all((candidate>=lo)&(candidate<=hi)):solutions.append(candidate)
            for k in range(2):
                j=1-k
                for value in [lo,hi]:
                    c=np.zeros(2);c[k]=value;c[j]=np.clip((rhs[j]-f[j,k]*value)/f[j,j],lo,hi);solutions.append(c)
            candidate=min(solutions,key=lambda c:float(c.dot(f).dot(c)-2*c.dot(rhs)))
        r=y-np.einsum('nlj,j->nl',b,candidate);radius=np.sqrt(np.sum(r*r,axis=1))
        newweight=base*np.minimum(1.,3*scale/np.maximum(radius,1e-20))
        done=np.max(abs(candidate-coef))<1e-9
        coef=candidate;weight=newweight
        if done:break
    f=np.einsum('n,nij->ij',weight,h);diag=np.sqrt(np.diag(f));eigen=np.linalg.eigvalsh(f/diag[:,None]/diag[None,:])
    return dict(coef=coef,iterations=iteration+1,normal_matrix=f,scaled_eigenvalues=eigen,
                downweighted=int(((weight<.999*base)&(base>0)).sum()),residual=r,weight=weight)

def main():
    with np.load(OUT/'measured.npz') as f:a={k:f[k] for k in ['run','event','R','fold','t','t_native']}
    with np.load(OUT/'moments.npz') as f:g={k:f[k] for k in ['run','event','valid','x','moment']}
    if not np.array_equal(a['run'],g['run']) or not np.array_equal(a['event'],g['event']):raise ValueError('Identity mismatch')
    if not g['valid'].all():raise ValueError('Invalid geometry')
    p,w=projector(g['x']);train=a['fold']!=0;test=~train
    train_runs,inverse=np.unique(a['run'][train],return_inverse=True);eval_runs=np.unique(a['run'][test])
    if np.intersect1d(train_runs,eval_runs).size:raise ValueError('Run leakage')
    masses=[12*.9314940954-6*.0005109989461,13.00335483507*.9314940954-6*.0005109989461]
    models=[];evaluations=[];bootstrap=[];rng=np.random.RandomState(10311);replicates=64
    boot_weights=[]
    for r in range(replicates):boot_weights.append(np.bincount(rng.randint(len(train_runs),size=len(train_runs)),minlength=len(train_runs))[inverse])
    for mass in masses:
        u=np.sqrt(1+(mass/(6*a['R']))**2);b=np.array([m.first_order(u[e],mass,g['moment'][e]) for e in range(len(u))])
        pb=np.einsum('nlk,nkj->nlj',p,b)
        for time_name in ['t_native','t']:
            t=a[time_name];y=np.einsum('nij,nj->ni',p,t)
            radius=np.linalg.norm(y[train],axis=1);scale=float(np.median(radius)/np.sqrt(2*np.log(2)))
            unbounded=fit(pb[train],y[train],scale)
            bounded=fit(pb[train],y[train],scale,bounds=(0.,3.))
            variants={'unbounded':unbounded,'nonnegative_box_0_3':bounded}
            tag=('native' if time_name=='t_native' else 'newcalib')+'_'+('C12_design' if mass==masses[0] else 'C13_design')
            entry=dict(tag=tag,mass_design_GeV=mass,robust_scale_ns=scale,variants={name:{k:v for k,v in f.items() if k not in ['residual','weight']} for name,f in variants.items()})
            freeze=OUT/(tag+'_contrast_freeze.json');freeze.write_text(json.dumps(native(entry),indent=2,allow_nan=False)+'\n')
            # Freeze is written before the run holdout is evaluated.
            result=dict(tag=tag,freeze_sha256=hashlib.sha256(freeze.read_bytes()).hexdigest(),variants={})
            for name,model in variants.items():
                coef=model['coef'];correction=np.einsum('nlj,j->nl',b[test],coef)
                residual=y[test]-np.einsum('nlk,nk->nl',p[test],correction)
                shift=-np.sum(w[test]*correction,axis=1)
                result['variants'][name]=dict(projected_time_baseline=[stats(y[test,l]) for l in range(4)],
                    projected_time_corrected=[stats(residual[:,l]) for l in range(4)],u_shift=stats(shift),
                    first_order_validity_warning=bool(np.max(abs(coef))>2),beta_resolution_measured=False)
            boot_coef=[]
            for bw in boot_weights:boot_coef.append(fit(pb[train],y[train],scale,base_weight=bw)['coef'])
            boot_coef=np.asarray(boot_coef);sensitivity=np.einsum('nl,nlj->nj',w[test],b[test]);shift_samples=-boot_coef.dot(sensitivity.T)
            boot=dict(tag=tag,replicates=replicates,coefficient_percentiles=np.percentile(boot_coef,[2.5,50,97.5],axis=0),
                u_calibration_sampling_sd=stats(np.std(shift_samples,axis=0,ddof=1)),coefficient_samples=boot_coef,
                interpretation='Whole calibration-run resampling; unbounded contrast diagnostic; fixed external timing, rigidity and material assumptions not resampled')
            models.append(entry);evaluations.append(result);bootstrap.append(boot)
            print('FIT',tag,'scale',scale,'unbounded',unbounded['coef'],'bounded',bounded['coef'],
                  'eigen',unbounded['scaled_eigenvalues'],'bootstrap',boot['coefficient_percentiles'],flush=True)
    report=dict(complete=True,status='diagnostic_only_not_a_deployable_calibration',train_events=int(train.sum()),train_runs=len(train_runs),
        evaluation_events=int(test.sum()),evaluation_runs=len(eval_runs),models=models,evaluation=evaluations,bootstrap=bootstrap,
        constraints=['Existing timing offsets frozen','Only two shared mean material scales','Both event t0 and inverse speed projected out',
            'No MC or RICH inputs','No final event isotope assignment','No old model overwritten'],
        omissions=['Rigidity error convolution and unknown isotope mixture','Timing calibration covariance','Full nonlinear mean loss in fitting',
            'Normalized selected-observation likelihood and detector noise model','Passive straggling, active deposit response and nuclear tails'],
        box_note='0<=kappa<=3 is a declared finite diagnostic envelope, not a proved physical range or optimum',
        RICH_read=False,MC_read=False,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'contrast_diagnostic.json').write_text(json.dumps(native(report),indent=2,allow_nan=False)+'\n')
    print('DIAGNOSTIC_COMPLETE',flush=True)

if __name__=='__main__':main()
