#!/usr/bin/env python3
"""Energy-only closure with an already calibrated clock baseline.

The generator has zero residual clock offsets by construction. No truth
clock correction is given to the fit. This isolates the simpler problem
from the earlier test that reopened weak clock directions.
"""
import json
import numpy as np
import flight_calibration_core as core
import research_flight_calibration_closure as toy
import research_flight_calibration_anchor as anchor


def generate(domain,n,seed):
    rng=np.random.RandomState(seed)
    x=np.array([-2.17,-2.04,2.04,2.17])[None,:]*rng.uniform(1.,1.3,(n,1))
    x+=rng.normal(0,.018,(n,4));x-=x.mean(axis=1)[:,None]
    epoch=rng.randint(0,toy.PERIODS,n)
    high=rng.uniform(size=n)<.3
    sensitivity=np.where(high,.001,1.)
    amplitude=np.clip(rng.normal(0,1,(n,4)),-2.,2.)
    position=np.stack([.5*(x[:,0]+x[:,1]),.5*x[:,1],.5*x[:,2],.5*(x[:,2]+x[:,3])],axis=1)
    b=np.stack([core.impulse_kernel(x,position[:,k])*(sensitivity*(1+.3*amplitude[:,k]))[:,None] for k in range(4)],axis=-1)
    u=np.where(high,rng.uniform(1.0002,1.003,n),rng.uniform(1.15,1.65,n))
    sigma=.045 if domain==0 else .065
    theta=toy.parameters(domain)[:4]
    t0=rng.normal(0,1,n)
    t=t0[:,None]+u[:,None]*x+np.einsum('nmp,p->nm',b,theta)+rng.normal(0,sigma,x.shape)
    return dict(x=x,times=t,design=b,period=epoch,high=high),dict(u=u,theta=theta,sigma=sigma)


def main():
    rows=[]
    # All seeds are reported; no successful seed is selected.
    for domain in [0,1]:
        for seed in [308,309,310,311]:
            m,_=generate(domain,50000,seed+domain*100)
            v,truth=generate(domain,12000,seed+2000+domain*100)
            fit=toy.train(m)
            path=toy.OUT/('fixed_timing_domain_%d_seed_%d_freeze.json'%(domain,seed))
            path.write_text(json.dumps(dict(theta=fit['theta'].tolist(),sigma=fit['estimated_sigma_ns'],
                core_sha256=toy.digest(core.__file__),script_sha256=toy.digest(__file__),
                clock_baseline='Zero residual clock offsets in the generator, held fixed in fit',truth_used=False),indent=2)+'\n')
            pred=toy.prediction(v,fit['theta'])
            raw=toy.prediction(v,np.zeros(4))
            known=toy.prediction(v,truth['theta'])
            spread=anchor.calibration_error_scale(v,fit['covariance_on_estimable_subspace']*fit['estimated_sigma_ns']**2)
            rows.append(dict(domain=domain,seed=seed,rank=fit['data_rank'],
                calibration_u_sigma_median=float(np.median(spread)),
                metrics={k:toy.stats(p-truth['u']) for k,p in [('raw',raw),('energy_only_fit',pred),('known_generator_parameters',known)]}))
    result=dict(complete=True,study='Control with frozen correct time baseline; standalone synthetic problem',
        RICH_used=False,flight_fit=False,MC_response_transferred=False,rows=rows,
        limitations='Correct linearized energy template and independent noise assumed. No unmodeled clock residual, energy-response drift, passive-loss fluctuations, or beta-dependent measurement error in covariates.',
        script_sha256=toy.digest(__file__),core_sha256=toy.digest(core.__file__))
    (toy.OUT/'fixed_timing_closure.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
