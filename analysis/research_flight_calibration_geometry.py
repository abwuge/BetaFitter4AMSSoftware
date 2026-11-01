#!/usr/bin/env python3
"""Geometry-only information audit on measured flight training events.

Impulse locations below are diagnostic templates, not a material map or a
calibration of energy losses. No RICH or MC field is opened.
"""
import json
import hashlib
from pathlib import Path
import numpy as np
import flight_calibration_core as core

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/flight_calibration_research'
ALLOW = {'run','event','source','u','R','relR','bin','continuous','bars'}


def qs(v):
    return np.percentile(v, [5,25,50,75,95]).tolist()


def eigen_summary(b, x):
    r = core.solve_global(x, np.zeros_like(x), b, 1., allow_rank_deficient=True)
    return dict(rank=r['data_rank'], parameter_count=b.shape[2],
                scaled_eigenvalues=r['scaled_eigenvalues'].tolist())


def main():
    path=ROOT/'results/time_local_beta/training.npz'
    with np.load(str(path), allow_pickle=False) as f:
        assert set(f.files)==ALLOW
        run=f['run']; r=f['R']; n=len(run)
        idx=np.random.RandomState(20261012).choice(n, min(n,100000), replace=False)
        cc=f['continuous'][idx].astype(float)
        r=r[idx]; run=run[idx]
    x=cc[:,10:14].copy()
    orient=np.sign(x[:,3]-x[:,0]); x*=orient[:,None]
    assert (np.diff(x,axis=1)>0).all()
    xx=np.stack([.5*(x[:,0]+x[:,1]),.5*x[:,1],.5*x[:,2],.5*(x[:,2]+x[:,3])],axis=1)
    h=np.stack([core.impulse_kernel(x,xx[:,k]) for k in range(4)],axis=-1)
    projected=core.project(x,h,1.)
    event_singular=np.linalg.svd(projected,compute_uv=False)
    central=core.impulse_kernel(x,0.)
    hp=core.project(x,central,1.)
    hc=central-central.mean(axis=1)[:,None]
    ratio=np.linalg.norm(hp,axis=1)/np.linalg.norm(hc,axis=1)
    # 50 ps per layer is a declared illustrative noise model, not an estimate.
    sigma=.05
    per_event_sigma_u_jump=sigma/np.linalg.norm(hp,axis=1)
    rows=[]
    for lo,hi in [(2.5,6),(6,20),(20,200)]:
        m=(r>=lo)&(r<hi)
        alpha=11.17793**2/(6*r[m])**3
        sensitivity_energy=alpha*np.linalg.norm(hp[m],axis=1)
        rows.append(dict(low=lo,high=hi,n=int(m.sum()),
            sigma_total_energy_GeV_for_common_parameter_illustration=float(sigma/np.sqrt(np.sum(sensitivity_energy**2))),
            median_single_event_energy_sigma_GeV=float(np.median(sigma/sensitivity_energy))))
    # Exact ambiguity of independent four-layer variances at fixed geometry.
    v=x[0]
    q=np.linalg.qr(np.column_stack([np.ones(4),v]),mode='complete')[0][:,2:].T
    t=np.stack([q[0]**2,q[1]**2,np.sqrt(2)*q[0]*q[1]])
    _,sv,vt=np.linalg.svd(t)
    null=vt[-1]/np.max(np.abs(vt[-1]))
    variance_a=np.ones(4)
    variance_b=variance_a+.8*null
    cov_a=(q*variance_a).dot(q.T); cov_b=(q*variance_b).dot(q.T)
    weights=(v-v.mean())/np.sum((v-v.mean())**2)
    covariance_ambiguity=dict(x_ns=v.tolist(),contrast_variance_map_rank=int(np.linalg.matrix_rank(t)),
        singular_values=sv.tolist(),variance_a=variance_a.tolist(),variance_b=variance_b.tolist(),
        max_contrast_covariance_difference=float(np.max(np.abs(cov_a-cov_b))),
        slope_variance_a=float(np.dot(weights**2,variance_a)),slope_variance_b=float(np.dot(weights**2,variance_b)),
        slope_contrast_covariance_a=q.dot(variance_a*weights).tolist(),
        slope_contrast_covariance_b=q.dot(variance_b*weights).tolist())
    # Distinct global effects may be identifiable because covariates vary.
    energy=cc[:,18:29]
    proxy=np.column_stack([np.mean(energy[:,:2],axis=1),np.mean(energy[:,4:8],axis=1),
                           np.mean(energy[:,8:],axis=1),np.mean(energy[:,2:4],axis=1)])
    proxy=(proxy-np.median(proxy,axis=0))/np.maximum(np.std(proxy,axis=0),1e-9)
    proxy=np.clip(proxy,-4,4)
    results=dict(complete=True,training_pool_n=n,sample_n=len(x),RICH_used=False,MC_used=False,
        source_training_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        path_ns_quantiles_by_layer=[qs(x[:,i]) for i in range(4)],
        station_gap_over_full_span_quantiles=[qs(np.diff(x,axis=1)[:,i]/np.ptp(x,axis=1)) for i in range(3)],
        event_rank_two_check_max_third_singular=float(np.max(event_singular[:,2])),
        impulse_amplitude_four_global_parameters=eigen_summary(h,x),
        impulse_times_measured_energy_proxies=eigen_summary(h*proxy[:,None,:],x),
        central_loss_projected_to_known_speed_norm_ratio=qs(ratio),
        single_event_inverse_speed_jump_sigma_at_50ps=qs(per_event_sigma_u_jump),
        illustrative_energy_information=rows,
        noise_assumption='Independent equal 0.05 ns layer noise only for illustrative sensitivity numbers. No measured noise calibration or CRLB for the full nonlinear experiment.',
        reference='Timing coordinate mean=0, not a claim of exact L5 reference. Projection is invariant to choice of reference origin.',
        template_assumption='Unit inverse-speed impulses halfway in S1-S2, S2-origin, origin-S3, S3-S4. Energy conversion uses nominal C12 mass and measured R only for information scaling.',
        proxy_limit='Existing beta-fed tk_ql plus TOF deposits are used solely to test algebraic variation of columns; rank is not proof of independent energy information.',
        covariance_ambiguity=covariance_ambiguity)
    OUT.mkdir(exist_ok=True)
    (OUT/'geometry_information.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps({k:results[k] for k in ['sample_n','station_gap_over_full_span_quantiles','impulse_amplitude_four_global_parameters','impulse_times_measured_energy_proxies','central_loss_projected_to_known_speed_norm_ratio','single_event_inverse_speed_jump_sigma_at_50ps','illustrative_energy_information','covariance_ambiguity']},indent=2),flush=True)

if __name__=='__main__': main()
