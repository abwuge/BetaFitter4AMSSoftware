#!/usr/bin/env python3
"""Freeze free-slope and continuous-mass transport outputs before truth scoring."""
import json
import subprocess
import sys
from pathlib import Path
import numpy as np
import flight_conditional_likelihood as f
import free_timing_transport as tr
import flight_mean_loss as phys
from research_conditional_mc import OUT,ROOT,inputs,stamp
from research_flight_conditional import write,digest,load

LABELS=['native_gaussianR','native_studentR5']


def main():
    if (OUT/'output_freeze.json').exists():raise FileExistsError('Output protocol is already frozen')
    a,g=inputs();test=a['fold']==0;x=g['x'][test];t=a['t_native'][test];R=a['R'][test];moment=g['moment'][test]
    check=subprocess.run([sys.executable,str(ROOT/'tests/test_free_timing_transport.py')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    (OUT/'transport_unit_tests.log').write_text(check.stdout)
    if check.returncode:raise RuntimeError(check.stdout)
    m0,m1=tr.compress(moment);xc=x-x.mean(1)[:,None];ols=xc/np.sum(xc*xc,axis=1)[:,None]
    s0=np.sum(ols*t,axis=1);pred={k:a[k][test] for k in ['run','event','source','entry','R','fold']}
    pred.update(ols_u=s0,constant_u=np.ones(len(x)),native_u=1/a['tof_betah'][test])
    ols_transport=tr.solve_center(s0,R,ols,moment,[1,1])
    if not ols_transport['converged'].all():raise ValueError('OLS transport inverse failed')
    pred['ols_unit_transport_u']=ols_transport['u']
    records=[]
    for label in LABELS:
        folder=OUT/label;fr=json.loads((folder/'freeze.json').read_text());ev=json.loads((folder/'measured_evaluation.json').read_text())
        if ev['freeze_sha256']!=digest(folder/'freeze.json'):raise ValueError('Changed model freeze')
        pp=load(folder/'predictions_before_truth.npz')
        if not np.array_equal(pp['run'],pred['run']) or not np.array_equal(pp['event'],pred['event']):raise ValueError('Prediction identity mismatch')
        for name,key in [('fitted','model'),('zero_loss','null_model'),('unit_material','unit_material_model')]:
            theta=np.asarray(fr[key]['theta']);w=tr.weights(theta,x);s=np.sum(w*t,axis=1)
            if not np.allclose(s,pp[name+'_u'],atol=1e-13,rtol=0):raise ValueError('Free-slope replay differs')
            tag=label+'_'+name;pred[tag+'_u']=s
            scales=theta[:2];center=tr.solve_center(s,R,w,moment,scales)
            if not center['converged'].all():raise ValueError('Transport inverse failed')
            pred[tag+'_transport_u']=center['u']
            # Forward perturbations include the physical material-time change.
            ucheck=np.linspace(.96,1.45,len(x));kernel=tr.continuous_kernel(ucheck,R,m0,m1)
            forward=ucheck[:,None]*x+np.einsum('nlg,g->nl',kernel,scales)+np.linspace(-1,1,len(x))[:,None]
            closure=tr.solve_center(np.sum(w*forward,axis=1),R,w,moment,scales)
            injection_error=float(np.max(abs(closure['u']-ucheck)))
            if injection_error>1e-10:raise ValueError('Physical forward/inverse closure failed')
            step=1e-4;up=tr.solve_center(s+step,R,w,moment,scales)['u'];dn=tr.solve_center(s-step,R,w,moment,scales)['u']
            derivative=(up-dn)/(2*step)
            records.append(dict(label=tag,kappa=scales,correction_percentiles=np.percentile(s-center['u'],[0,16,50,84,100]),
                max_inverse_residual=float(np.max(abs(center['residual']))),physical_injection_max_error=injection_error,
                stored_time_slope_response_percentiles=np.percentile(derivative,[0,16,50,84,100]),
                beta_above_one_before=int(((s>0)&(s<1)).sum()),beta_above_one_after=int(((center['u']>0)&(center['u']<1)).sum()),
                mass_rule='Continuous m=6*R*sqrt(u_L5^2-1), no isotope posterior',freeze_sha256=digest(folder/'freeze.json')))
    np.savez_compressed(OUT/'outputs_before_truth.npz',**pred)
    # This protocol is set before any reference branch is opened by this code.
    protocol=dict(created_utc=stamp(),complete=True,events=len(x),runs=len(np.unique(pred['run'])),
        records=records,output_sha256=digest(OUT/'outputs_before_truth.npz'),
        model_freeze_sha256={label:digest(OUT/label/'freeze.json') for label in LABELS},
        module_sha256=digest(Path(tr.__file__)),test_sha256=digest(ROOT/'tests/test_free_timing_transport.py'),
        truth_read=False,RICH_read=False,MC_coefficients_transferred_to_flight=False,
        reference_rules=dict(primary_plane='Native MC acceptance plane10; same generated particle, trackID, charge and mass; positive finite momentum and near-L5 coordinate',
            native_association='Run/event and stored acceptance values checked against raw; primary reference within 1 cm of measured L5 track state',
            reference_development_status='Existing MC development sources, with measured-only fitting newly frozen; not new simulation or blind flight data',
            unsupported_raw='Stored acceptance-reference evaluation is reported separately when raw MC cannot be audited'),
        metrics=['inverse-beta error bias, median, W68=(q84-q16)/2, RMSE, absolute-tail >0.02 and >0.05',
                 'beta error bias, median, W68, RMSE, absolute-tail >0.02 and >0.05',
                 'paired entire-evaluation-run resampling with frozen outputs; training-run refit uncertainty separately'],
        bins=dict(measured_R=[2.5,3,3.5,4.5,6],reference_beta=[0,.8,.85,.9,.95,1.01]),
        parameter_selection='No candidate is selected or retuned using the reference metrics',
        transport_scope='First-order material timing inverse, continuous event mass, observed R at z=0 approximated as L5 momentum; beta>1 extension retains measured output',
        limitations=['Conditional model selection truncation and source-dependent noise remain assumptions',
                     'Mean transport uses existing material map and elemental Bethe expression, not a fitted visible-energy response',
                     'Free output response is exact for affine time shifts; material-corrected output is tested against full physical forward time changes'])
    write(OUT/'output_freeze.json',protocol)
    print('MC_OUTPUTS_FROZEN',len(x),'events',len(records),'branches',flush=True)
    for row in records:print(row['label'],'kappa',row['kappa'],'du percentiles',row['correction_percentiles'],'slope response',row['stored_time_slope_response_percentiles'],flush=True)


if __name__=='__main__':main()
