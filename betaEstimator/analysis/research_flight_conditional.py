#!/usr/bin/env python3
"""Numerical checks and flight-only fits of the conditional mean-loss model."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import flight_conditional_likelihood as f
import flight_mean_loss as physical

ROOT=f.ROOT;BASE=ROOT/'results/flight_calibration_research/mean_loss_likelihood';OUT=f.OUT

def load(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k] for k in z.files}

def inputs():
    with np.load(BASE/'measured.npz',allow_pickle=False) as z:a={k:z[k] for k in ['run','event','R','relR','fold','t_native','t']}
    g=load(BASE/'moments.npz')
    if not np.array_equal(a['run'],g['run']) or not np.array_equal(a['event'],g['event']):raise ValueError('Identity mismatch')
    return a,g

def native(x):
    if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [native(v) for v in x]
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    return x

def write(p,x):Path(p).write_text(json.dumps(native(x),indent=2,allow_nan=False)+'\n')
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def density_integral(model,theta,index,nodes):
    info=model.numpy_event(theta,index);means=info['means'].reshape(-1,3);weights=info['conditional_rigidity_weights'].ravel()
    center=weights.dot(means);difference=means-center;l,p=f.precision(theta);sigma=np.linalg.inv(p)
    cov=5/3*sigma+np.einsum('n,ni,nj->ij',weights,difference,difference);s=np.linalg.cholesky(cov)
    x,w=np.polynomial.legendre.leggauss(nodes);z=np.tan(np.pi*.5*x);jw=w*(np.pi*.5)*(1+z*z)
    grid=np.stack(np.meshgrid(z,z,z,indexing='ij'),axis=-1).reshape(-1,3)
    volume=np.einsum('i,j,k->ijk',jw,jw,jw).ravel()*np.linalg.det(s)
    points=center+grid.dot(s.T);total=0.
    normal=float(f.gammaln(4)-f.gammaln(2.5)-1.5*np.log(5*np.pi)+np.log(np.diag(l)).sum())
    for start in range(0,len(points),256):
        xx=points[start:start+256];delta=xx[:,None,:]-means[None,:,:]
        radius=np.sum(np.matmul(delta,l)**2,axis=-1)
        pdf=np.exp(normal-4*np.log1p(radius/5)).dot(weights)
        total+=float(np.dot(pdf,volume[start:start+len(xx)]))
    return total

def verify():
    a,g=inputs();ix=np.flatnonzero(a['fold']!=0)[::300][:32];theta=f.initial(.7,.25)
    theta[1]=1.3;theta[5:8]=[.17,-.11,.13];theta[8:10]=[1.2,-.3];theta[11]=.25;theta[12]=.1
    checks=[]
    for rn in [0.,5.]:
        model=f.ConditionalModel(a,g,ix,nodes=192,rigidity_nu=rn)
        loss,gradient,events=model.evaluate(theta,events=True)
        reference=np.array([-model.numpy_event(theta,e)['loglike'] for e in range(len(ix))])
        errors=[]
        for k in range(13):
            step=1e-5;up=theta.copy();down=theta.copy();up[k]+=step;down[k]-=step
            finite=(model.evaluate(up)[0]-model.evaluate(down)[0])/(2*step)
            errors.append(float(finite-gradient[k]))
        if np.max(abs(reference-events))>1e-11 or np.max(np.abs(errors))>2e-7:raise ValueError('Density or gradient verification failure')
        checks.append(dict(rigidity_nu=rn,numpy_cpp_max_nll_error=float(np.max(abs(reference-events))),gradient_errors=errors,max_gradient_error=float(np.max(np.abs(errors)))))
    model=f.ConditionalModel(a,g,ix,nodes=192)
    kernel_errors=[]
    for species,mass in enumerate(f.MASSES):
        for qindex in [5,80,140]:
            direct=f.D.dot(physical.first_order(model.u[species,qindex],mass,g['moment'][ix[0]]))
            compressed=model.moment0[0]*model.a[species,qindex]+model.moment1[0]*model.b[species,qindex]
            kernel_errors.append(float(np.max(abs(direct-compressed))))
    if max(kernel_errors)>1e-12:raise ValueError('Material compression failed')
    integral={str(n):density_integral(model,theta,0,n) for n in [20,32,48]}
    if abs(integral['48']-1)>2e-5:raise ValueError('Conditional density integration failed: '+str(integral))
    normal=sum(model.numpy_event(theta,0)['conditional_rigidity_weights'].ravel())
    record=dict(passed=True,gradient_checks=checks,material_kernel_max_error=max(kernel_errors),density_integrals=integral,
        conditional_latent_weights_sum=float(normal),dimension=3,density_measure='Orthonormal time contrasts in ns cubed',
        scope='Numerical normalization and gradients under the stated conditional model; not detector-model validation',
        source_sha256={name:digest(ROOT/'analysis'/name) for name in ['flight_conditional_likelihood.py','flight_conditional_likelihood.cpp','research_flight_conditional.py']},
        library_sha256=digest(OUT/'libconditional.so'))
    write(OUT/'numerical_verification.json',record);print(json.dumps(record,indent=2),flush=True)

def optimize(model,start,fixed=None,label='fit',maxiter=600):
    fixed={} if fixed is None else fixed;free=np.array([k for k in range(13) if k not in fixed]);initial=np.asarray(start,float).copy()
    for k,v in fixed.items():initial[k]=v
    calls=[0];t0=time.time()
    def expand(x):
        out=initial.copy();out[free]=x;return out
    def objective(x):
        theta=expand(x);loss,grad=model.evaluate(theta);calls[0]+=1
        if calls[0]%100==0:
            print(label,'eval',calls[0],'NLL',loss,'kappa',theta[:2],'seconds',round(time.time()-t0,1),flush=True)
        return loss,grad[free]
    result=minimize(objective,initial[free],jac=True,method='L-BFGS-B',bounds=[f.BOUNDS[k] for k in free],
                    options=dict(maxiter=maxiter,maxfun=3000,ftol=2e-11,gtol=2e-7,maxls=40))
    theta=expand(result.x);loss,grad=model.evaluate(theta)
    pg=grad.copy()
    for k,(lo,hi) in enumerate(f.BOUNDS):
        if k in fixed or (theta[k]<=lo+1e-7 and grad[k]>0) or (theta[k]>=hi-1e-7 and grad[k]<0):pg[k]=0
    return dict(theta=theta,nll=loss,success=bool(result.success),message=str(result.message),iterations=int(result.nit),evaluations=int(result.nfev),
        elapsed_seconds=time.time()-t0,gradient=grad,projected_gradient_max=float(np.max(abs(pg))),
        fixed={f.PARAMETERS[k]:v for k,v in fixed.items()},boundaries=[f.PARAMETERS[k] for k,(lo,hi) in enumerate(f.BOUNDS) if abs(theta[k]-lo)<1e-5 or abs(theta[k]-hi)<1e-5])

def fit_model(time_field,rigidity_nu,starts,nodes):
    check=json.loads((OUT/'numerical_verification.json').read_text())
    if not check['passed']:raise ValueError('Numerical checks not passed')
    a,g=inputs();train=a['fold']!=0;test=~train
    if np.intersect1d(a['run'][train],a['run'][test]).size:raise ValueError('Run overlap')
    label=('native' if time_field=='t_native' else 'newcalib')+('_gaussianR' if rigidity_nu<=0 else '_studentR%d'%rigidity_nu)
    folder=OUT/label;folder.mkdir(exist_ok=True)
    trainmodel=f.ConditionalModel(a,g,train,time_field=time_field,nodes=nodes,rigidity_nu=rigidity_nu)
    protocol=dict(events=int(train.sum()),runs=len(np.unique(a['run'][train])),evaluation_events=int(test.sum()),evaluation_runs=len(np.unique(a['run'][test])),
        time_field=time_field,rigidity_nu=rigidity_nu,time_student_nu=5,quadrature_nodes=nodes,true_inverse_rigidity_range=trainmodel.qrange,
        parameters=f.PARAMETERS,bounds=f.BOUNDS,masses_GeV=f.MASSES.tolist(),RICH_read=False,MC_read=False,
        density='p(Dt | measured inverse rigidity, stored inverse-rigidity error, fitted geometry, selected sample)',
        physics='First-order mean material timing; two isotope masses integrated; both spectrum shape and isotope odds fitted',
        limitations=['Selected-sample conditional noise assumption, not full readout-selection likelihood','Two-coefficient log-spectrum and linear log-isotope-odds shared across sources',
            'First-order rather than full nonlinear mean transport in fitting','No active deposit response or passive straggling',
            'Fixed native timing reference with no external calibration covariance','Finite latent-rigidity domain and parameter bounds require sensitivity checks'],
        inference='This is a calibration likelihood only; its event isotope/spectrum posterior is not a deployed beta estimator')
    write(folder/'protocol.json',protocol)
    fits=[]
    choices=[(1.,.1),(0.,.5),(2.,.9)][:starts]
    for i,(kappa,fraction) in enumerate(choices):
        result=optimize(trainmodel,f.initial(kappa,fraction),label=label+'_start%d'%i)
        fits.append(result);write(folder/('start_%d.json'%i),result)
        print('START_DONE',label,i,'nll',result['nll'],'theta',result['theta'],'gradient',result['projected_gradient_max'],flush=True)
    best=min(fits,key=lambda r:r['nll']);best['selection_rule']='Minimum training NLL among all completed prescribed starts'
    null=optimize(trainmodel,best['theta'],fixed={0:0.,1:0.},label=label+'_null')
    fixed_one=optimize(trainmodel,best['theta'],fixed={0:1.,1:1.},label=label+'_fixed1')
    freeze=dict(model=best,null_model=null,unit_material_model=fixed_one,starts=fits,protocol_sha256=digest(folder/'protocol.json'),
        module_sha256=digest(ROOT/'analysis/flight_conditional_likelihood.py'),cpp_sha256=digest(ROOT/'analysis/flight_conditional_likelihood.cpp'),library_sha256=digest(OUT/'libconditional.so'))
    write(folder/'freeze.json',freeze)
    evaluation={};predictions={k:a[k][test] for k in ['run','event']}
    for nq in [nodes,2*nodes]:
        model=f.ConditionalModel(a,g,test,time_field=time_field,nodes=nq,rigidity_nu=rigidity_nu)
        row={}
        for name,fit in [('fitted',best),('zero_loss',null),('unit_material',fixed_one)]:
            nll,gradient,events=model.evaluate(fit['theta'],events=True);row[name]=nll
            predictions['nll_%s_%d'%(name,nq)]=events
        evaluation[str(nq)]=row
    np.savez_compressed(folder/'evaluation_nll.npz',**predictions)
    write(folder/'evaluation.json',dict(complete=True,events=int(test.sum()),runs=len(np.unique(a['run'][test])),mean_conditional_nll=evaluation,
        freeze_sha256=digest(folder/'freeze.json'),metric_scope='Held-run conditional timing density, not beta-resolution measurement'))
    print('FIT_COMPLETE',label,evaluation,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['verify','fit']);p.add_argument('--time-field',choices=['t_native','t'],default='t_native');p.add_argument('--rigidity-nu',type=float,default=0.);p.add_argument('--starts',type=int,default=3);p.add_argument('--nodes',type=int,default=128);args=p.parse_args()
    if args.stage=='verify':verify()
    else:fit_model(args.time_field,args.rigidity_nu,args.starts,args.nodes)
