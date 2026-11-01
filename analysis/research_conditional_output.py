#!/usr/bin/env python3
"""Check held-run likelihood uncertainty and free-slope timing-output sensitivity.

No truth or RICH measurements enter this diagnostic. Reported u differences
are differences between reconstructions, not their errors or resolutions.
"""
import json
from pathlib import Path
import numpy as np
import flight_conditional_likelihood as f
from research_flight_conditional import inputs,write,digest

OUT=f.OUT
LABELS=['native_gaussianR','newcalib_gaussianR','native_studentR5']

def stats(a):
    a=np.asarray(a,float).ravel();q=np.percentile(a,[16,50,84])
    return dict(n=len(a),mean=float(a.mean()),median=float(q[1]),w68=float((q[2]-q[0])/2),rms=float(np.sqrt(np.mean(a*a))),max_abs=float(np.max(abs(a))))

def weights(theta,x):
    _,precision=f.precision(theta);dx=np.asarray(x,float).dot(f.D.T)
    pa=dx.dot(precision);den=np.sum(pa*dx,axis=1)
    return pa.dot(f.D)/den[:,None]

def hessian(model,theta,step):
    columns=[]
    for k in range(13):
        eps=step*max(1.,abs(theta[k]));up=theta.copy();dn=theta.copy();up[k]+=eps;dn[k]-=eps
        columns.append((model.evaluate(up)[1]-model.evaluate(dn)[1])/(2*eps))
    raw=np.asarray(columns).T
    return (raw+raw.T)*.5,float(np.max(abs(raw-raw.T)))

def main():
    a,g=inputs();train=a['fold']!=0;test=~train;x=g['x'][test]
    run=a['run'][test];runs,inverse=np.unique(run,return_inverse=True);counts=np.bincount(inverse)
    rng=np.random.RandomState(12102026);boot=np.array([np.bincount(rng.randint(len(runs),size=len(runs)),minlength=len(runs)) for _ in range(512)])
    den=boot.dot(counts);base_x=x-x.mean(1)[:,None];ols=base_x/np.sum(base_x*base_x,axis=1)[:,None]
    records=[];saved_weights={};predictions=dict(run=run,event=a['event'][test])
    for label in LABELS:
        folder=OUT/label;protocol=json.loads((folder/'protocol.json').read_text());freeze=json.loads((folder/'freeze.json').read_text())
        evaluation=json.loads((folder/'evaluation.json').read_text())
        if evaluation['freeze_sha256']!=digest(folder/'freeze.json'):raise ValueError('Freeze changed')
        theta=np.asarray(freeze['model']['theta'],float);time_field=protocol['time_field'];t=a[time_field][test]
        with np.load(folder/'evaluation_nll.npz',allow_pickle=False) as nll:
            if not np.array_equal(nll['run'],run) or not np.array_equal(nll['event'],a['event'][test]):raise ValueError('Evaluation identity mismatch')
            diffs={name:nll['nll_fitted_256']-nll['nll_'+name+'_256'] for name in ['zero_loss','unit_material']}
        likelihood={}
        for name,delta in diffs.items():
            sums=np.bincount(inverse,weights=delta);samples=boot.dot(sums)/den
            likelihood[name]=dict(mean_nll_difference=float(delta.mean()),total_nll_difference=float(delta.sum()),
                percentile_2p5_50_97p5=np.percentile(samples,[2.5,50,97.5]).tolist(),
                convention='Fitted minus nuisance-refitted comparison; negative favors fitted density',
                uncertainty='512 paired held-run resamples with all calibration parameters frozen; no retraining or model uncertainty')
        variants={};wfit=weights(theta,x);saved_weights[label]=wfit
        for name,key in [('fitted','model'),('zero_loss','null_model'),('unit_material','unit_material_model')]:
            th=np.asarray(freeze[key]['theta'],float);w=weights(th,x);u=np.sum(w*t,axis=1)
            predicted_u=np.sum(ols*t,axis=1)
            injection_a=np.linspace(-2,2,len(t));injection_b=np.linspace(-.2,.2,len(t))
            injected=np.sum(w*(t+injection_a[:,None]+injection_b[:,None]*x),axis=1)-u-injection_b
            variants[name]=dict(median_weights_ns_inverse=np.median(w,axis=0).tolist(),
                max_constant_response_error=float(np.max(abs(w.sum(axis=1)))),
                max_slope_response_error=float(np.max(abs(np.sum(w*x,axis=1)-1))),
                max_injected_response_error=float(np.max(abs(injected))),
                u_difference_from_OLS=stats(u-predicted_u),
                u_difference_from_fitted_weights=stats(np.sum((w-wfit)*t,axis=1)),
                measured_beta_above_one=int(((u>0)&(u<1)).sum()),nonpositive_inverse_speed=int((u<=0).sum()),
                estimator='Free time slope with frozen global precision; no event isotope posterior, spectrum prior, or mean-loss correction')
            predictions[label+'_'+name+'_u']=u
            if np.max(abs(injected))>1e-10:raise ValueError('Free-slope response failure')
        # A velocity-parallel covariance addition is an example of a harmless
        # noise parameter direction. This is not a proof for every weak mode.
        _,p=f.precision(theta);sigma=np.linalg.inv(p);gauge_errors=[]
        for e in range(min(64,len(x))):
            dx=f.D.dot(x[e]);added=sigma+.03**2*np.outer(dx,dx)
            inv=np.linalg.inv(added);wnew=f.D.T.dot(inv).dot(dx)/dx.dot(inv).dot(dx)
            gauge_errors.append(np.max(abs(wnew-wfit[e])))
        model=f.ConditionalModel(a,g,train,time_field=time_field,nodes=128,rigidity_nu=protocol['rigidity_nu'])
        h,asym=hessian(model,theta,2e-4);fine,asymfine=hessian(model,theta,1e-4)
        nn=fine[2:,2:];eig=np.linalg.eigvalsh(nn)
        eff=fine[:2,:2]-fine[:2,2:].dot(np.linalg.solve(nn,fine[2:,:2]))
        # An interior quadratic diagnostic is useful for scale, but boundary
        # minima and misspecified response models forbid calling it a CI.
        hinfo=dict(mean_loss_nuisance_profiled_hessian_total=(len(model.y)*eff).tolist(),
            total_profiled_eigenvalues=np.linalg.eigvalsh(len(model.y)*eff).tolist(),
            nuisance_mean_hessian_eigenvalues=eig.tolist(),
            maximum_hessian_step_difference=float(np.max(abs(h-fine))),maximum_antisymmetric_component=asymfine,
            boundaries=freeze['model']['boundaries'],
            interpretation='Observed local smooth-extension curvature at a bounded fit, not Fisher information, sampling covariance or confidence intervals')
        if eig.min()>0 and np.linalg.eigvalsh(eff).min()>0:
            hinfo['interior_quadratic_scale_only_kappa_sd']=np.sqrt(np.diag(np.linalg.inv(len(model.y)*eff))).tolist()
        records.append(dict(label=label,events=int(test.sum()),runs=len(runs),kappa=theta[:2].tolist(),
            likelihood=likelihood,free_slope_output=variants,local_curvature=hinfo,
            parallel_covariance_gauge_max_weight_difference=float(max(gauge_errors)),freeze_sha256=digest(folder/'freeze.json')))
        print('OUTPUT',label,'kappa',theta[:2],'NLLvs0',likelihood['zero_loss'],'curvature',hinfo,flush=True)
    cross=np.sum((saved_weights['native_studentR5']-saved_weights['native_gaussianR'])*a['t_native'][test],axis=1)
    result=dict(complete=True,records=records,paired_resamples=512,
        native_rigidity_noise_family_weight_output_difference=stats(cross),
        limitations=['No beta truth or independent flight reference used','Free-slope output is not yet corrected to physical L5 speed',
            'Coherent injection applies to the stored time interface, not upstream beta-fed reconstruction',
            'Hessian diagnostics at parameter bounds do not supply calibrated confidence intervals',
            'Same fitted covariance across all source epochs, charge patterns and geometries remains a model assumption',
            'Parallel covariance invariance is one algebraic example; other weak directions may alter beta'],
        RICH_read=False,MC_read=False,script_sha256=digest(__file__),module_sha256=digest(Path(f.__file__)))
    write(OUT/'output_diagnostic.json',result);np.savez_compressed(OUT/'free_slope_predictions.npz',**predictions)
    print('OUTPUT_COMPLETE','native_noise_family_u_difference',stats(cross),flush=True)

if __name__=='__main__':main()
