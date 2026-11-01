#!/usr/bin/env python3
"""Joint-bar sensitivity and measured-response profiles; no beta refit."""
from datetime import datetime,timezone
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import attribute_tofcalib_time as t
from prepare_tofcalib_validation import OUT


def main():
    t.w.check_freeze()
    with np.load(str(OUT/'predictions.npz')) as f:t.P={k:f[k] for k in f.files}
    with np.load(str(OUT/'diagnostic_covariates.npz')) as f:t.C={k:f[k] for k in f.files}
    for k in ('run','event','source'):assert np.array_equal(t.P[k],t.C[k])
    with np.load(str(t.w.DEST/'predictions.npz')) as f:old=f['run'][f['reference_mask']]
    t.REF=t.P['reference_mask']&(t.P['run']>=old.min())&(t.P['run']<=old.max())
    years=np.arange(2011,2027);limits=[datetime(int(y),1,1,tzinfo=timezone.utc).timestamp() for y in years]
    t.YEAR=years[np.searchsorted(limits,t.P['run'],side='right')-1];t.DIMS=t.build_dimensions()
    bars=np.sum(t.C['bars'].astype(np.int64)*np.array([1,10,100,1000]),axis=1)
    rb=np.searchsorted(t.w.R_EDGES,t.P['R'],side='right')-1
    joint=(bars*len(t.w.R_EDGES)+rb)*2+t.P['NaF'].astype(int)
    _,code=np.unique(joint,return_inverse=True)
    additional=t.DIMS[3]+[('joint_four_bars_R_radiator',code)]
    e0=1/t.P['raw_beta_rich']-1/t.P['rich_beta'];e1=1/t.P['new_beta_rich']-1/t.P['rich_beta']
    delta=t.P['raw_u']-t.P['corrected_u']
    assert np.max(abs((e0-delta)-e1))<1e-12
    ref0=np.flatnonzero(t.REF);rows=[]
    def moments(x,v):
        mean=float(np.average(x,weights=v));var=float(np.average((x-mean)**2,weights=v));return mean,var
    for year in np.unique(t.YEAR):
        target0=np.flatnonzero((t.YEAR==year)&t.P['reference_mask'])
        ref,target=t.supported(ref0,target0,t.DIMS[3]);weights,balance=t.balance(ref,target,t.DIMS[3])
        m0,v0=moments(e0[target],weights);m1,v1=moments(e1[target],weights);md,vd=moments(delta[target],weights)
        cov=float(np.average((e0[target]-m0)*(delta[target]-md),weights=weights))
        assert np.isclose(v1-v0,vd-2*cov,atol=1e-14,rtol=1e-10)
        low,high=t.quantile(e0[target],weights,[.05,.95]);central=(e0[target]>=low)&(e0[target]<=high)
        e=e0[target][central];d=delta[target][central];cw=weights[central]
        me,ve=moments(e,cw);mc,vc=moments(d,cw);co=float(np.average((e-me)*(d-mc),weights=cw))
        row=dict(year=int(year),n=len(target),kinematics_quality_geometry_balanced=t.measure(target,weights),
            balance=balance,mean_log_edep=[float(np.average(t.C['continuous'][target,i],weights=weights)) for i in range(18,22)],
            mean_log_tracker=[float(np.average(t.C['continuous'][target,i],weights=weights)) for i in range(22,29)],
            timing_input_medians=[float(t.quantile(t.C['continuous'][target,i],weights,[.5])[0]) for i in range(14,18)],
            inverse_beta_variance_decomposition=dict(raw_variance=v0,new_variance=v1,correction_variance=vd,
                twice_error_correction_covariance=2*cov,closure_error=(v1-v0)-(vd-2*cov)),
            central_90pct_correlation_diagnostic=dict(n=int(central.sum()),error_correction_correlation=co/np.sqrt(ve*vc),
                correction_sd=float(np.sqrt(vc)),raw_error_sd=float(np.sqrt(ve)),variance_cancellation_ratio=2*co/vc))
        if int(year) in [2011,2014,2015,2021,2023,2025]:
            ar,at=t.supported(ref0,target0,additional);aw,ab=t.balance(ar,at,additional)
            row['joint_bar_sensitivity']=dict(reference=t.measure(ar,np.ones(len(ar))),target=t.measure(at,aw),
                balance=ab,reference_coverage=len(ar)/len(ref0),target_coverage=len(at)/len(target0))
        rows.append(row)
        print(year,'core correlation',round(row['central_90pct_correlation_diagnostic']['error_correction_correlation'],3),
              'jointbars',row.get('joint_bar_sensitivity',{}).get('target',{}).get('gain_pct'),flush=True)
    t.w.ev.dump(OUT/'response_attribution.json',dict(complete=True,beta_model_refitted=False,rows=rows,
        core_diagnostic_selection='Central 90% of baseline inverse-beta residual, for covariance diagnosis only; main gains always use complete matched populations',
        causal_limit='Calendar association and measured response, not identification of a specific calibration release or physical aging mechanism',
        predictions_sha256=t.w.ev.digest(OUT/'predictions.npz'),covariates_sha256=t.w.ev.digest(OUT/'diagnostic_covariates.npz'),
        script_sha256=t.w.ev.digest(Path(__file__)),balancing_script_sha256=t.w.ev.digest(Path(t.__file__))))
    fig,axes=plt.subplots(2,2,figsize=(11,8));years=[r['year'] for r in rows]
    for key,label in [('old_w68','Linear TOF'),('new_w68','Frozen correction')]:
        axes[0,0].plot(years,[1000*r['kinematics_quality_geometry_balanced'][key] for r in rows],'o-',label=label)
    axes[0,0].set_ylabel('Beta residual W68 x 1000');axes[0,0].legend(frameon=False)
    for j in range(4):axes[0,1].plot(years,[r['mean_log_edep'][j] for r in rows],'o-',label='Layer %d'%(j+1))
    axes[0,1].set_ylabel('Mean log(TOF energy deposit / 36)');axes[0,1].legend(frameon=False,fontsize=8)
    axes[1,0].plot(years,[r['central_90pct_correlation_diagnostic']['error_correction_correlation'] for r in rows],'o-')
    axes[1,0].set_ylabel('Core inverse-beta error / correction correlation')
    for j in range(4):axes[1,1].plot(years,[r['timing_input_medians'][j] for r in rows],'o-',label='Layer %d'%(j+1))
    axes[1,1].set_ylabel('Median TOF time residual / time span');axes[1,1].legend(frameon=False,fontsize=8)
    for ax in axes.flat:ax.set_xlabel('Calendar year (UTC)');ax.set_xticks(np.arange(2011,2026,3))
    fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(str(OUT/'figures'/('response_attribution.'+ext)),dpi=170)
    plt.close(fig)


if __name__=='__main__':main()
