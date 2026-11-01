#!/usr/bin/env python3
"""Verify mean-loss kernels against bidirectional cell integration, without labels."""
import json
import hashlib
from pathlib import Path
import numpy as np
import uproot
import flight_mean_loss as m

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/flight_calibration_research/mean_loss_likelihood'

def stats(a):
    a=np.asarray(a,float).ravel()
    if not np.isfinite(a).all():raise ValueError('Nonfinite propagation result')
    return dict(n=len(a),mean=float(a.mean()),median=float(np.median(a)),rms=float(np.sqrt(np.mean(a*a))),max_abs=float(np.max(abs(a))),abs_p95=float(np.percentile(abs(a),95)))

def main():
    with np.load(OUT/'measured.npz') as f:a={k:f[k] for k in ['run','event','R']}
    with np.load(OUT/'moments.npz') as f:g={k:f[k] for k in ['x','moment','valid']}
    with uproot.open(OUT/'geometry.root') as f:raw=f['geometry'].arrays(library='np')
    # Uniform ranks in measured rigidity cover the full declared range.
    order=np.argsort(a['R']);ix=order[np.unique(np.linspace(0,len(order)-1,128).astype(int))]
    # Atomic masses from NIST; legacy CODATA unit/electron values are recorded,
    # not fitted. Neglecting atomic binding is an explicit test approximation.
    mass_unit_GeV=.9314940954;electron_GeV=.0005109989461
    masses=[12*mass_unit_GeV-6*electron_GeV,13.00335483507*mass_unit_GeV-6*electron_GeV]
    configs=[(1.,1.),(2.,2.),(1.,0.),(0.,1.)]
    result=[];zero=[];linear_derivative=[];first_order_slope=[];nonlinear_slope=[];run=[];event=[];reference_R=[]
    for mass in masses:
        for scales in configs:
            time_errors=[];slope_errors=[];convergence=[];losses=[];first_slope=[];full_slope=[]
            for e in ix:
                hi=np.asarray(raw['high'][e],float);lo=np.asarray(raw['low'][e],float);le=np.asarray(raw['length'][e],float)
                col=np.asarray(raw['abundance'][e],float).reshape(-1,9)*m.ATOMIC_MASS
                group=np.isin(raw['segment'][e],[1,2]).astype(np.int64)
                args=(hi,lo,le,col,group,np.asarray(raw['at_high'][e],np.int64),np.asarray(raw['at_low'][e],np.int64),float(raw['center_z'][e]),raw['tof_z'][e].astype(float))
                u=np.sqrt(1+(mass/(6*a['R'][e]))**2);x=g['x'][e];z=x-x.mean();w=z/np.dot(z,z)
                kappa=np.asarray(scales,float)
                full,loss=m.propagate_cells(u,mass,kappa,*args,substeps=1)
                finer,_=m.propagate_cells(u,mass,kappa,*args,substeps=4)
                b=m.first_order(u,mass,g['moment'][e]);linear=x*u+b.dot(kappa)
                time_errors.extend(linear-finer);slope_errors.append(np.dot(w,linear-finer));convergence.append(np.dot(w,full-finer));losses.append(loss)
                first_slope.append(np.dot(w,b.dot(kappa)));full_slope.append(np.dot(w,finer-x*u))
                if scales==(1.,1.):
                    no_loss,_=m.propagate_cells(u,mass,np.zeros(2),*args);zero.extend(no_loss-x*u)
                    tiny,_=m.propagate_cells(u,mass,np.full(2,1e-4),*args,substeps=4)
                    linear_derivative.extend((tiny-x*u)/1e-4-b.sum(1))
                    run.append(a['run'][e]);event.append(a['event'][e]);reference_R.append(a['R'][e]);first_order_slope.append(first_slope[-1]);nonlinear_slope.append(full_slope[-1])
            result.append(dict(mass_GeV=mass,kappa=scales,events=len(ix),first_order_minus_nonlinear_time_ns=stats(time_errors),
                first_order_minus_nonlinear_u=stats(slope_errors),one_vs_four_substep_u=stats(convergence),
                first_order_speed_correction=stats(first_slope),nonlinear_speed_correction=stats(full_slope),
                integrated_loss_mean_GeV=np.mean(losses,axis=0).tolist(),integrated_loss_max_GeV=np.max(losses,axis=0).tolist()))
    report=dict(complete=True,geometry_events=len(a['R']),tested_events=len(ix),tests_per_event=len(masses)*len(configs),
        zero_loss_time_ns=stats(zero),small_loss_derivative_time_ns=stats(linear_derivative),cases=result,
        mass_definition=dict(atomic_relative_masses=[12.,13.00335483507],atomic_mass_unit_GeV=mass_unit_GeV,electron_mass_GeV=electron_GeV,
            atomic_binding_omitted=True,source='NIST Atomic Weights and Isotopic Compositions for Carbon; unit conversion from installed scipy legacy CODATA values'),
        scope='Numerical consistency under the same Bethe/material model, not validation of actual detector stopping power or resolution',
        reference='Measured rigidity used only to choose numerical test states; not asserted to equal the event true momentum',
        RICH_read=False,MC_read=False,fit_performed=False)
    if report['zero_loss_time_ns']['max_abs']>1e-10:raise ValueError('Zero-loss limit failed')
    if report['small_loss_derivative_time_ns']['max_abs']>1e-4:raise ValueError('First-order derivative mismatch')
    report['source_sha256']={name:hashlib.sha256((ROOT/'analysis'/name).read_bytes()).hexdigest() for name in ['flight_mean_loss.py','verify_flight_mean_loss.py']}
    (OUT/'propagation_verification.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(OUT/'propagation_check.npz',run=run,event=event,R=reference_R,first_order_slope=first_order_slope,nonlinear_slope=nonlinear_slope)
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
