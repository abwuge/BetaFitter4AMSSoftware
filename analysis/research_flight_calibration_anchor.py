#!/usr/bin/env python3
"""Test a global high-speed calibration anchor, absent at inference.

The high-speed mean of 1 is a deliberately approximate physics anchor in
this toy. Its bias is measured, not hidden. This is not a flight fit.
"""
import json
import hashlib
from pathlib import Path
import numpy as np
import flight_calibration_core as core
import research_flight_calibration_closure as toy


def fit_anchored(m,sigma,mean,anchor_sigma):
    x,t,b=m['x'],m['times'],m['design']
    w=1/sigma**2
    a=np.stack([np.ones_like(x),x],axis=2)
    local=w*np.einsum('nmi,nmj->nij',a,a)
    lam=m['high'].astype(float)/anchor_sigma**2
    local[:,1,1]+=lam
    local_inv=np.linalg.inv(local)
    cross=w*np.einsum('nmi,nmp->nip',a,b)
    rhs_local=w*np.einsum('nmi,nm->ni',a,t)
    rhs_local[:,1]+=lam*mean
    fisher=w*np.einsum('nmp,nmq->pq',b,b)-np.einsum('nip,nij,njq->pq',cross,local_inv,cross)
    rhs=w*np.einsum('nmp,nm->p',b,t)-np.einsum('nip,nij,nj->p',cross,local_inv,rhs_local)
    covariance=np.linalg.inv(fisher)
    theta=covariance.dot(rhs)
    scales=np.sqrt(np.diag(fisher))
    eigen=np.linalg.eigvalsh(fisher/scales[:,None]/scales[None,:])
    return dict(theta=theta,covariance=covariance,scaled_eigenvalues=eigen)


def calibration_error_scale(m,covariance):
    zero=np.zeros(m['design'].shape[2])
    base=core.reconstruct(m['x'],m['times'],m['design'],1.,zero)
    jac=-np.einsum('nm,nmp->np',base['slope_weights'],m['design'])
    v=np.einsum('np,pq,nq->n',jac,covariance,jac)
    return np.sqrt(np.maximum(v,0))


def main():
    out=toy.OUT
    results=[]
    for domain in [0,1]:
        m,_=toy.generate(domain,toy.PROTOCOL['train_events_per_domain'],toy.PROTOCOL['seed']+1000+domain)
        v,truth=toy.generate(domain,toy.PROTOCOL['test_events_per_domain'],toy.PROTOCOL['seed']+2000+domain)
        naive=toy.train(m)
        sigma=naive['estimated_sigma_ns']
        anchored=fit_anchored(m,sigma,mean=1.,anchor_sigma=.003)
        settings=dict(high_speed_mean=1.,high_speed_spread=.003,
                      uses_per_event_truth=False,anchor_used_at_inference=False,
                      anchor_bias='Generator high-speed mean is above 1; this intentionally approximate anchor can shift the global beta scale.')
        frozen=dict(theta=anchored['theta'].tolist(),settings=settings,
                    estimated_sigma_ns=sigma,core_sha256=toy.digest(core.__file__),
                    script_sha256=toy.digest(__file__))
        (out/('anchored_domain_%d_freeze.json'%domain)).write_text(json.dumps(frozen,indent=2)+'\n')
        raw=toy.prediction(v,np.zeros(toy.P))
        unanchored=toy.prediction(v,naive['theta'])
        pred=toy.prediction(v,anchored['theta'])
        sn=calibration_error_scale(v,naive['covariance_on_estimable_subspace']*sigma**2)
        sa=calibration_error_scale(v,anchored['covariance'])
        shifted=dict(v,times=v['times']+.01*v['x']-3.)
        record=dict(domain=domain,settings=settings,
            data_only_scaled_eigenvalues=naive['scaled_eigenvalues'].tolist(),
            anchored_scaled_eigenvalues=anchored['scaled_eigenvalues'].tolist(),
            calibration_induced_u_sigma_median_unanchored=float(np.median(sn)),
            calibration_induced_u_sigma_median_anchored=float(np.median(sa)),
            metrics={k:toy.stats(p-truth['u']) for k,p in [('raw',raw),('unanchored',unanchored),('anchored',pred)]},
            affine_response_error=float(np.max(np.abs(toy.prediction(shifted,anchored['theta'])-pred-.01))))
        results.append(record)
    result=dict(complete=True,study='Posthoc repair of the failed pure-contrast toy; no flight/RICH model selection',
        rows=results,script_sha256=toy.digest(__file__),toy_sha256=toy.digest(toy.__file__),
        interpretation='A physical ensemble anchor can stabilize GLOBAL calibration while inference retains free event slopes. Approximate anchor bias remains. This is not an unbiasedness or optimality proof for real data.')
    (out/'anchored_closure.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
