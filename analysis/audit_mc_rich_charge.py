#!/usr/bin/env python3
"""Independent C12 check of RICH charge quality using native primary brackets."""
import argparse
import json
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import uproot
import awkward as ak
import periodic_joint_adc as period
import rich_transport as material
from audit_low_rich_residuals import summary
from research_flight_conditional import ROOT,load,write,digest
from research_conditional_mc import stamp

POOL=period.OUT.parent/'additional_mc'
OUT=period.OUT/'mc_rich_charge_audit'
FIELDS=['run','event','itrtrack','irich','rich_itrtrack','rich_beta','rich_ebeta1','rich_pb','rich_used','rich_q','rich_npe',
    'rich_good','rich_clean','rich_NaF','rich_pmt','rich_PMTChargeConsistency','rich_BetaConsistency','rich_pos','rich_tile']


def plan():
    OUT.mkdir(exist_ok=True)
    if (OUT/'protocol.json').exists():raise FileExistsError('MC ring-reference audit already declared')
    write(OUT/'protocol.json',dict(created_utc=stamp(),source_sha256=digest(Path(__file__)),
        measured_sha256=digest(POOL/'measured.npz'),source_manifest_sha256=digest(POOL/'protocol.json'),
        native_manifest_sha256=digest(POOL/'reference/complete.json'),quality_source_sha256=digest(Path(material.__file__)),
        fields=FIELDS,quality='Unchanged strict rich_transport quality, NaF, plus original measured Tracker total chi-square in (0,10].',
        charge_windows=[[4.5,7.5],[5.,7.],[5.5,6.5]],
        primary='Exact generator particle, track, parent, mass and charge at both stored native planes 17 and 18, bracketing the actual stored RICH emission z; unique primary records and unchanged measured Tracker/BetaH association.',
        geometric_acceptance='Both primary-state straight extrapolations to stored RICH z within 1 cm of its fitted intersection; downward directions; finite positive momentum and carbon charge.',
        speed_reference='Use the two primary momenta as an explicit interval, not a fabricated exact momentum at the RICH plane. Report bracket widths and only bound-robust incompatibilities. Monotone-loss interpretation assumes no energy-gaining interaction between the surfaces.',
        midpoint_use='On interval widths <=0.001 only, midpoint diagnostics have a deterministic half-width error bound; not labeled exact truth.',
        incompatibility_distance=.01,narrow_bracket_max=.001,
        exposure='Existing additional C12 MC measured pool and native primary records; independent of flight calibration. No C13 generalization claimed.',
        no_reference_to_reconstruction_training=True,no_MC_coefficients_to_flight=True,no_truth_error_selection=True,
        purpose='Validate the mechanism and possible bias of reference-only charge windows before any new flight confirmatory use. Report every window, including failures.'))
    print('MC_RICH_CHARGE_PLAN',flush=True)


def read():
    p=json.loads((OUT/'protocol.json').read_text())
    if digest(Path(__file__))!=p['source_sha256'] or digest(Path(material.__file__))!=p['quality_source_sha256']:
        raise ValueError('MC ring audit implementation changed')
    if (OUT/'readout.npz').exists():raise FileExistsError('MC RICH fields already read')
    a=load(POOL/'measured.npz');n=len(a['run']);records=json.loads((POOL/'protocol.json').read_text())['records']
    complete=json.loads((POOL/'reference/complete.json').read_text())
    if digest(POOL/'measured.npz')!=p['measured_sha256'] or digest(POOL/'reference/complete.json')!=p['native_manifest_sha256']:
        raise ValueError('MC primary inputs changed')
    def one(record):
        source=record['source'];ix=np.flatnonzero(a['source']==source);entry=a['entry'][ix]
        with uproot.open(record['path']) as stream:
            raw=stream['amstreea'].arrays(FIELDS,entry_stop=record['entry_stop'],library='np')
        rich={key:value[entry] for key,value in raw.items()}
        for key in ['run','event','itrtrack']:
            if not np.array_equal(rich[key],a[key][ix]):raise ValueError('MC ring event association differs')
        path=POOL/'reference/parts'/('%05d.root'%source)
        if digest(path)!=complete['part_hashes'][str(source)]:raise ValueError('Native primary reference changed')
        with uproot.open(path) as stream:native=ak.to_list(stream['reference'].arrays(library='ak'))
        lookup={(row['run'],row['event']):row for row in native}
        planes=np.full((len(ix),2,8),np.nan);reason=np.full(len(ix),'missing_primary',dtype='U40')
        for index,i in enumerate(ix):
            row=lookup[(int(a['run'][i]),int(a['event'][i]))]
            if not (row['found'] and row['track_exists'] and row['beta_exists'] and row['beta_track']==row['requested_track']):
                reason[index]='native_association_failure';continue
            if row['requested_track']!=a['itrtrack'][i] or row['requested_beta']!=a['ibetah'][i]:raise ValueError('Primary target differs')
            chosen=[]
            for plane in [17,18]:
                candidates=[j for j,label in enumerate(row['plane']) if label==plane and row['particle'][j]==row['gen_particle'] and row['track'][j]==row['gen_track'] and row['parent'][j]==row['gen_parent'] and row['mass'][j]==row['gen_mass'] and row['charge'][j]==row['gen_charge']]
                if len(candidates)!=1:break
                chosen.append(candidates[0])
            if len(chosen)!=2:continue
            reason[index]='records_found'
            for k,j in enumerate(chosen):
                planes[index,k]=np.r_[row['momentum'][j],row['mass'][j],row['position'][3*j:3*j+3],row['direction'][3*j:3*j+3]]
            if row['gen_charge']!=6:reason[index]='primary_charge_not_carbon'
        return ix,rich,planes,reason,dict(source=source,path=record['path'],events=len(ix),native_part_sha256=digest(path))
    output={};native=np.full((n,2,8),np.nan);reason=np.full(n,'unread',dtype='U40');seen=np.zeros(n,bool);evidence=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for count,(ix,rich,planes,why,note) in enumerate(pool.map(one,records),1):
            for key,value in rich.items():
                if key not in output:output[key]=np.empty((n,)+value.shape[1:],value.dtype)
                output[key][ix]=value
            native[ix]=planes;reason[ix]=why;seen[ix]=True;evidence.append(note)
            if count%32==0:print('MC_RICH_CHARGE_READ',count,len(records),flush=True)
    if not seen.all():raise ValueError('Incomplete MC ring audit')
    output.update(primary_bracket=native,primary_reason=reason,source=a['source'],entry=a['entry'])
    np.savez_compressed(OUT/'readout.npz',**output)
    write(OUT/'readout_manifest.json',dict(complete=True,created_utc=stamp(),events=n,records=evidence,
        protocol_sha256=digest(OUT/'protocol.json'),readout_sha256=digest(OUT/'readout.npz')))


def evaluate():
    p=json.loads((OUT/'protocol.json').read_text());manifest=json.loads((OUT/'readout_manifest.json').read_text())
    if (OUT/'evaluation.json').exists():raise FileExistsError('MC charge-reference score already saved')
    amendment=json.loads((OUT/'scoring_amendment.json').read_text())
    if amendment['original_source_sha256']!=p['source_sha256'] or digest(Path(__file__))!=amendment['corrected_source_sha256']:
        raise ValueError('MC scoring key correction differs from its amendment')
    if digest(OUT/'protocol.json')!=amendment['protocol_sha256'] or digest(OUT/'readout.npz')!=manifest['readout_sha256']:
        raise ValueError('MC charge audit input changed')
    a=load(POOL/'measured.npz');r=load(OUT/'readout.npz');v=r['primary_bracket'];n=len(a['run'])
    chi=a['chi'][:,2];basic=material.rich_quality(r,strict=True)&r['rich_NaF']&np.isfinite(chi)&(chi>0)&(chi<=10)
    z=r['rich_pos'][:,2];position=v[:,:,2:5];direction=v[:,:,5:8];momentum=v[:,:,0];mass=v[:,:,1]
    with np.errstate(divide='ignore',invalid='ignore'):
        extrapolated=position[:,:,:2]+(z[:,None]-position[:,:,2])[:,:,None]*direction[:,:,:2]/direction[:,:,2:3]
        distance=np.linalg.norm(extrapolated-r['rich_pos'][:,None,:2],axis=2)
        beta=momentum/np.hypot(momentum,mass)
    bracket=(r['primary_reason']=='records_found')&np.isfinite(v).all((1,2))&(momentum>0).all(1)&(mass>0).all(1)
    bracket&=(direction[:,:,2]<-.05).all(1)&(position[:,0,2]>z)&(position[:,1,2]<z)&(distance<=1.).all(1)
    bracket&=(momentum[:,0]>=momentum[:,1])&(mass[:,0]==mass[:,1])
    lower=beta[:,1];upper=beta[:,0];width=upper-lower;midpoint=.5*(upper+lower);rich_beta=r['rich_beta'][:,0]
    discrepancy=np.maximum(lower-rich_beta,np.maximum(rich_beta-upper,0))
    narrow=bracket&(width<=p['narrow_bracket_max'])
    charge=np.sqrt(np.maximum(r['rich_q'][:,0],0));valid_charge=np.isfinite(r['rich_q'][:,0])&(r['rich_q'][:,0]>=0)
    cohorts=[('original_quality',basic)]+[('charge_%g_%g'%(lo,hi),basic&valid_charge&(charge>=lo)&(charge<=hi)) for lo,hi in p['charge_windows']]
    rows=[]
    for name,mask in cohorts:
        known=mask&bracket;thin=mask&narrow;count=int(known.sum());bad=known&(discrepancy>p['incompatibility_distance'])
        row=dict(cohort=name,events=int(mask.sum()),primary_bracket_valid=count,narrow_bracket_events=int(thin.sum()),
            bracket_width=summary(width[known]),robust_incompatible=int(bad.sum()),
            robust_incompatible_fraction=float(bad.sum()/count) if count else None,
            midpoint_residual_narrow_only=summary((rich_beta-midpoint)[thin]),
            midpoint_pointwise_error_bound=summary(.5*width[thin]),
            charge=summary(charge[known]),bins=[])
        for lo,hi in [(0,.78),(.78,.82),(.82,.86),(.86,.9),(.9,.94),(.94,1.01)]:
            bx=known&(midpoint>=lo)&(midpoint<hi);bn=bx&narrow
            row['bins'].append(dict(primary_beta_low=lo,primary_beta_high=hi,events=int(bx.sum()),
                robust_incompatible=int((bad&bx).sum()),narrow_events=int(bn.sum()),midpoint_residual=summary((rich_beta-midpoint)[bn])))
        rows.append(row)
        print('MC_RICH_CHARGE_CLOSURE',name,'quality',row['events'],'bracket',count,'incompatible',row['robust_incompatible'],
            'midpoint_bias',row['midpoint_residual_narrow_only'].get('mean'),flush=True)
    np.savez_compressed(OUT/'reference_diagnostic.npz',run=a['run'],event=a['event'],original_quality=basic,primary_valid=bracket,
        primary_lower_beta=lower,primary_upper_beta=upper,primary_distance=distance,rich_beta=rich_beta,rich_charge=charge,discrepancy_to_interval=discrepancy)
    write(OUT/'evaluation.json',dict(complete=True,created_utc=stamp(),events=n,primary_status_counts=dict(Counter(r['primary_reason'].tolist())),cohorts=rows,
        protocol_sha256=digest(OUT/'protocol.json'),readout_manifest_sha256=digest(OUT/'readout_manifest.json'),
        diagnostic_sha256=digest(OUT/'reference_diagnostic.npz'),scoring_amendment_sha256=digest(OUT/'scoring_amendment.json'),no_flight_models_updated=True,not_exact_RICH_plane_truth=True,
        limitation='C12 MC only; pre-existing mean-material and ring-response differences between simulation and flight are not calibrated by this check.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['plan','read','evaluate']);args=parser.parse_args();globals()[args.stage]()
