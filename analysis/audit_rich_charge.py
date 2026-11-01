#!/usr/bin/env python3
"""Read native-stored ring charge on the already exposed reference sources."""
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import uproot
from audit_low_rich_residuals import POOL,summary
import periodic_joint_adc as period
from research_flight_conditional import ROOT,load,write,digest
from research_conditional_mc import stamp
OUT=period.OUT/'rich_charge_audit'


def main():
    OUT.mkdir(exist_ok=True)
    if (OUT/'protocol.json').exists():raise FileExistsError('RICH charge audit already started')
    producer=ROOT.parent/'amsd69n/analysis_amsd69n.C'
    write(OUT/'protocol.json',dict(created_utc=stamp(),source_sha256=digest(Path(__file__)),producer_sha256=digest(producer),
        measured_sha256=digest(POOL/'measured.npz'),parent_sha256=digest(POOL/'reference/parent.npz'),
        reference_sha256=digest(POOL/'reference/reference.npz'),preparation_sha256=digest(POOL/'preparation.json'),
        field_definitions={'rich_q[0]':'ring->getCharge2Estimate()','rich_q[1]':'sqrt(ring->NpCol/ring->NpExp), or -1 for nonpositive expected count'},
        source_scope='Only parent sources already opened for the existing RICH evaluation.',
        hypothesis='Rings whose speed strongly disagrees with the upstream measured event may also have an incompatible ring charge estimate. This is a reference-consistency diagnostic, not truth adjudication.',
        diagnostic_groups=['all valid','R2.5-3GV all','R2.5-3GV final beta difference<-0.03','R2.5-3GV absolute beta difference<0.02'],
        outcome_grouping='Diagnostic labels only, not fitting inputs, not an event-rejection rule.',
        fixed_charge_windows=[[4.5,7.5],[5.,7.],[5.5,6.5]],
        windows_scope='Report every declared window without selecting one by beta improvement. No replacement reference or fresh-blind claim.',
        no_models_updated=True,no_production_selection_changed=True))
    a=load(POOL/'measured.npz');parent=load(POOL/'reference/parent.npz');ref=load(POOL/'reference/reference.npz')
    pred=load(period.OUT/'development_beta/flight/predictions_before_scoring.npz')
    prep=json.loads((POOL/'preparation.json').read_text());mapping={row['source']:row['path'] for row in prep['records']}
    def one(source):
        ix=np.flatnonzero(a['source']==source);entry=a['entry'][ix];low=int(entry.min());high=int(entry.max())+1
        with uproot.open(mapping[int(source)]) as f:
            raw=f['amstreea'].arrays(['run','event','itrtrack','irich','rich_q','rich_npe'],entry_start=low,entry_stop=high,library='np')
        selected={key:value[entry-low] for key,value in raw.items()}
        for key in ['run','event','itrtrack']:
            if not np.array_equal(selected[key],a[key][ix]):raise ValueError('Ring charge event association differs')
        for key in ['irich','rich_npe']:
            if not np.array_equal(selected[key],parent[key][ix]):raise ValueError('Ring charge parent reference differs')
        return ix,selected['rich_q'],dict(source=int(source),path=mapping[int(source)],events=len(ix),range=[low,high])
    charge=np.full((len(a['run']),2),np.nan);seen=np.zeros(len(charge),bool);sources=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for ix,value,note in pool.map(one,np.unique(a['source'])):
            charge[ix]=value;seen[ix]=True;sources.append(note)
    if not seen.all():raise ValueError('Ring charge read was incomplete')
    np.savez_compressed(OUT/'charge.npz',run=a['run'],event=a['event'],source=a['source'],entry=a['entry'],rich_q=charge)
    sqrt_charge=np.sqrt(np.maximum(charge[:,0],0));ratio=np.sqrt(np.maximum(parent['rich_npe'][:,0],0)/np.maximum(parent['rich_npe'][:,1],1e-10))
    qvalid=np.isfinite(charge).all(1)&(charge[:,0]>=0)&(charge[:,1]>=0)
    baseline=1/pred['newcalib_gaussianR_static_ols_fit']-ref['beta'];valid=ref['valid'];low=valid&(a['R']<3.)
    groups={'all_valid':valid,'low_R':low,'low_R_negative_discrepancy':low&(baseline<-.03),'low_R_core_diagnostic':low&(abs(baseline)<.02)}
    result=dict(complete=True,created_utc=stamp(),read_events=len(charge),source_files=sources,groups={},charge_windows=[],
        protocol_sha256=digest(OUT/'protocol.json'),charge_sha256=digest(OUT/'charge.npz'),no_events_removed_from_existing_evaluations=True)
    for name,mask in groups.items():
        result['groups'][name]=dict(events=int(mask.sum()),charge_valid=int((mask&qvalid).sum()),
            sqrt_getCharge2=summary(sqrt_charge[mask]),raw_np_charge=summary(charge[mask,1]),getter_npe_proxy=summary(ratio[mask]),
            getter_minus_raw=summary(sqrt_charge[mask]-charge[mask,1]),proxy_minus_raw=summary(ratio[mask]-charge[mask,1]))
        print('RICH_CHARGE_DIAGNOSTIC',name,int(mask.sum()),'median charge',result['groups'][name]['sqrt_getCharge2']['percentiles'][3],flush=True)
    for lower,upper in json.loads((OUT/'protocol.json').read_text())['fixed_charge_windows']:
        accepted=qvalid&(sqrt_charge>=lower)&(sqrt_charge<=upper)
        result['charge_windows'].append(dict(low=lower,high=upper,
            group_counts={name:dict(total=int(mask.sum()),inside=int((mask&accepted).sum())) for name,mask in groups.items()}))
    write(OUT/'evaluation.json',result)
    print('RICH_CHARGE_WINDOWS',json.dumps(result['charge_windows']),flush=True)


if __name__=='__main__':main()
