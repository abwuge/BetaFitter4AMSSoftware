#!/usr/bin/env python3
"""Diagnose, never remove, the low-rigidity flight-reference discrepancies."""
import json
from pathlib import Path
import numpy as np
import free_timing_transport as center
import periodic_joint_adc as period
from research_flight_conditional import ROOT,load,write,digest
from research_conditional_mc import stamp

POOL=period.OUT.parent/'prospective_flight'
OUT=period.OUT/'low_rigidity_reference_audit'


def summary(value):
    v=np.asarray(value,float);v=v[np.isfinite(v)]
    if not len(v):return dict(events=0)
    q=np.percentile(v,[0,2.5,16,50,84,97.5,100])
    return dict(events=len(v),mean=float(v.mean()),rms=float(np.sqrt(np.mean(v*v))),w68=float(.5*(q[4]-q[2])),percentiles=q)


def main():
    OUT.mkdir(exist_ok=True)
    if (OUT/'protocol.json').exists():raise FileExistsError('Low-rigidity diagnostic already started')
    sources={str(path.relative_to(ROOT)):digest(path) for path in [POOL/'measured.npz',POOL/'moments.npz',
        POOL/'reference/reference.npz',POOL/'reference/parent.npz',period.OUT/'development_beta/flight/predictions_before_scoring.npz']}
    write(OUT/'protocol.json',dict(created_utc=stamp(),source_sha256=digest(Path(__file__)),inputs=sources,
        purpose='Decompose the previously observed low-R beta residuals into measured slope, mean material propagation and reference variants. This is an exposed-reference diagnostic, not a new event-selection rule.',
        rigidity_range=[2.5,3.],reference_status='Previously opened reference pool; no new blind claim.',
        diagnostic_groups=['all reference-valid','R in [2.5,3)','R in [2.5,3) and final OLS beta-reference<-0.03','R in [2.5,3) and absolute final OLS residual<0.02'],
        residual_grouping='Used only to characterize failures. No fit, model choice, reported acceptance cohort or production cut may be changed by these diagnostic labels.',
        constant_mass_design='Use C12/C13 kinematic design only as a diagnostic; no event isotope assignment.',
        affected_fits=False,representative_events='Deterministic first ten run/event identities from the discrepant group, not the most extreme errors.'))
    a=load(POOL/'measured.npz');g=load(POOL/'moments.npz');ref=load(POOL/'reference/reference.npz');parent=load(POOL/'reference/parent.npz')
    pred=load(period.OUT/'development_beta/flight/predictions_before_scoring.npz')
    for source in [g,ref,parent,pred]:
        for key in ['run','event']:
            if not np.array_equal(source[key],a[key]):raise ValueError('Diagnostic event identities differ')
    n=len(a['run']);x=g['x'];xc=x-x.mean(1)[:,None];w=xc/np.sum(xc*xc,axis=1)[:,None]
    native_slope=np.sum(w*a['t_native'],axis=1);new_slope=np.sum(w*a['t'],axis=1)
    material=center.solve_center(new_slope,a['R'],w,g['moment'],np.array([2.,2.]))
    if not material['converged'].all():raise ValueError('L5 inverse failed in diagnostic')
    final=1/pred['newcalib_gaussianR_static_ols_fit'];rich=ref['beta'];error=final-rich
    truth_design={name:1/np.sqrt(1+(mass/(6*a['R']))**2) for name,mass in [('C12_design',11.174862),('C13_design',12.109483)]}
    quantities=dict(R=a['R'],relR=a['relR'],R_L5_minus_R=a['R_L5']-a['R'],stored_native_beta=a['tof_betah'],
        native_OLS_beta=1/native_slope,newcalib_OLS_beta=1/new_slope,L5_beta=1/material['u'],final_RICH_plane_beta=final,RICH_beta=rich,
        raw_native_minus_RICH=a['tof_betah']-rich,native_OLS_minus_RICH=1/native_slope-rich,
        newcalib_minus_RICH=1/new_slope-rich,L5_minus_RICH=1/material['u']-rich,final_minus_RICH=error,
        time_calibration_beta_change=1/new_slope-1/native_slope,L5_correction_beta=1/material['u']-1/new_slope,
        forward_RICH_beta_change=final-1/material['u'],RICH_error=ref['error'],
        RICH_Beta_minus_getBeta=ref['beta_variants'][:,1]-rich,RICH_BetaRefit_minus_getBeta=ref['beta_variants'][:,2]-rich,
        ring_probability=parent['rich_pb'],PMTs=parent['rich_pmt'],used_hits=parent['rich_used'],
        ring_consistency=parent['rich_BetaConsistency'],PMT_charge_consistency=parent['rich_PMTChargeConsistency'],
        collected_npe=parent['rich_npe'][:,0],expected_npe=parent['rich_npe'][:,1],
        ring_to_event_npe=parent['rich_npe'][:,0]/np.maximum(parent['rich_npe'][:,2],1e-10),
        rich_x=parent['rich_pos'][:,0],rich_y=parent['rich_pos'][:,1],native_time_chi=a['native_chi'][:,0],
        native_tracker_total_chi=a['native_chi'][:,2],charge=a['tk_qin'][:,0,2])
    quantities.update(truth_design)
    for key,value in truth_design.items():
        quantities['RICH_minus_'+key]=rich-value;quantities['native_minus_'+key]=a['tof_betah']-value
    valid=ref['valid'];low=valid&(a['R']>=2.5)&(a['R']<3.)
    groups={'all_reference_valid':valid,'low_R':low,'low_R_negative_discrepancy':low&(error<-.03),
        'low_R_core_diagnostic':low&(abs(error)<.02)}
    records={}
    for name,mask in groups.items():
        records[name]=dict(events=int(mask.sum()),runs=len(np.unique(a['run'][mask])),summary={key:summary(value[mask]) for key,value in quantities.items()})
    ix=np.flatnonzero(groups['low_R_negative_discrepancy']);order=np.lexsort((a['event'][ix],a['run'][ix]));ix=ix[order[:10]]
    examples=[]
    for i in ix:
        examples.append(dict(run=int(a['run'][i]),event=int(a['event'][i]),source=int(a['source'][i]),entry=int(a['entry'][i]),
            itrtrack=int(a['itrtrack'][i]),ibetah=int(a['ibetah'][i]),irich=int(parent['irich'][i]),
            values={key:float(value[i]) for key,value in quantities.items()}))
    write(OUT/'evaluation.json',dict(complete=True,created_utc=stamp(),groups=records,examples=examples,
        protocol_sha256=digest(OUT/'protocol.json'),no_events_removed=True,no_parameters_updated=True))
    for name,row in records.items():
        print('LOW_R_REFERENCE_AUDIT',name,row['events'],json.dumps({key:row['summary'][key] for key in ['native_OLS_minus_RICH','time_calibration_beta_change','L5_correction_beta','forward_RICH_beta_change','RICH_BetaRefit_minus_getBeta','native_minus_C12_design','RICH_minus_C12_design']},default=lambda a:a.tolist()),flush=True)


if __name__=='__main__':main()
