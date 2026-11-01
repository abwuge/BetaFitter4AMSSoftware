#!/usr/bin/env python3
"""Compare native event associations, timing values and available pulse sets."""
import hashlib
import json
from pathlib import Path
import numpy as np
import uproot

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/flight_calibration_research/native_timing_audit'

def native(v):
    if isinstance(v,dict):return {str(k):native(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [native(x) for x in v]
    if isinstance(v,np.ndarray):return native(v.tolist())
    if isinstance(v,np.generic):return native(v.item())
    if isinstance(v,float) and not np.isfinite(v):return None
    return v

def stats(v):
    v=np.asarray(v,float).ravel();v=v[np.isfinite(v)]
    if not len(v):return dict(n=0)
    return dict(n=len(v),median=float(np.median(v)),max_abs=float(np.max(abs(v))),rms=float(np.sqrt(np.mean(v*v))))

def pulses(a,e,l,s,q):
    o=int(a['pulse_offset'][e,l,s,q]);n=int(a['pulse_count'][e,l,s,q])
    return np.asarray(a['pulse_time'][e][o:o+n],float)

def main():
    protocol=json.loads((OUT/'protocol.json').read_text());rows=[]
    for record in protocol['records']:
        d=OUT/str(record['source'])
        with np.load(d/'original.npz') as f:b={k:f[k] for k in f.files}
        with uproot.open(d/'native.root') as f:a=f['timing'].arrays(library='np')
        identity={k:bool(np.array_equal(a[k],b[k])) for k in ['run','event','itrtrack','ibetah']}
        if not all(identity.values()) or not a['selected'].all() or not a['found'].all():raise ValueError('Native association failure')
        good=(a['goodtime']!=0)
        comparisons={
            'stored_time_vs_native_cluster':stats(a['time'][good]-b['tof_tl'][good]),
            'stored_length_vs_native_beta_length':stats(a['pathlen']-b['tof_leng']),
            'native_beta_time_vs_cluster_time':stats((a['beta_time']-a['time'])[good]),
            'stored_chi2_vs_native':stats(a['chi2']-b['tof_chist']),
            'used_count_matches_flags':bool(np.array_equal(a['usehit'],a['used'].sum(axis=1))),
            'stored_bar_matches_native':bool(np.array_equal(a['bar'],b['tof_barid']))}
        for q,key in enumerate(['tof_rtime','tof_rtime_ht','tof_rtime_sht']):comparisons[key]=stats(a['raw'][:,:,:,q]-b[key])
        ids=np.flatnonzero(a['event']==record['tail_event'])
        if len(ids)!=1:raise ValueError('Extreme identity not unique')
        e=int(ids[0]);extreme={k:a[k][e] for k in ['run','event','used','goodtime','pattern','status','sidebit','time','beta_time','etime','beta_etime','pathlen','residual','beta_value','beta_error','chi2','pulse_code','pulse_count','raw','corrected']}
        extreme['pulses_ns']=[[[pulses(a,e,l,s,q) for q in range(3)] for s in range(2)] for l in range(4)]
        extreme['pulse_time_convention']='FT-threshold time after native cable/slewing correction; signs are opposite particle arrival times; diagnostic only, no candidate chosen'
        controls=a['event']!=record['tail_event']
        row=dict(source=record['source'],events=len(a['run']),identity=identity,comparisons=comparisons,extreme=extreme,
            controls=dict(events=int(controls.sum()),four_used=int((a['usehit'][controls]==4).sum()),chi2=stats(a['chi2'][controls]),
                          multiple_LT_sides=int((a['pulse_count'][controls,:,:,0]>1).sum())),
            native_sha256=hashlib.sha256((d/'native.root').read_bytes()).hexdigest())
        rows.append(row)
    result=dict(complete=True,events=sum(r['events'] for r in rows),association_mismatches=0,rows=rows,
        RICH_read=False,MC_read=False,selection='Known extreme events plus same-run uniform controls; not an efficiency measurement',
        no_timing_correction_applied=True,training_membership=json.loads((OUT/'training_membership.json').read_text()))
    (OUT/'verification.json').write_text(json.dumps(native(result),indent=2,allow_nan=False)+'\n')
    for r in rows:
        print('SOURCE',r['source'],'COMPARISONS',r['comparisons'],flush=True)
        e=r['extreme'];print('EXTREME',native({k:e[k] for k in ['run','event','used','goodtime','time','beta_time','etime','beta_value','chi2','pulse_count']}),flush=True)
        print('LT_PULSES',native([[pulses_from[0] for pulses_from in layer] for layer in e['pulses_ns']]),flush=True)
    print('VERIFIED',result['events'],flush=True)

if __name__=='__main__':main()
