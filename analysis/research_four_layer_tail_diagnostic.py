#!/usr/bin/env python3
"""Retain explicit tail identities; do not turn residual size into a new cut."""
import json
import numpy as np
import uproot
import research_flight_timing_redundancy as r
BASE=r.ROOT/'results/flight_calibration_research'
with np.load(str(BASE/'four_layer_timing_native/measured.npz'),allow_pickle=False) as f:
    a={k:f[k] for k in f.files}
p,w=r.m.projector(a['x']);res=np.einsum('nij,nj->ni',p,a['t']);norm=np.linalg.norm(res,axis=1)
ix=np.flatnonzero(norm>10.)
records=json.loads((r.ROOT/'results/tofcalib_validation/preparation.json').read_text())['records']
lookup={int(v['source']):v['path'] for v in records};rows=[]
fields=['run','event','itrtrack','ibetah','tof_hsumh','tof_hsumhu','tof_chist','tof_chist_n','tof_barid','tof_tl','tof_tl_NewCalib','tof_rtime','tof_rtime_ht','tof_rtime_sht']
for i in ix:
    with uproot.open(lookup[int(a['source'][i])]) as f:
        tree=f['amstreea'];available=[key for key in fields if key in tree.keys()];b=tree.arrays(available,library='np')
    take=(b['run']==a['run'][i])&(b['event']==a['event'][i])
    if take.sum()!=1:raise ValueError('Tail identity mismatch')
    values={k:v[take][0] for k,v in b.items()}
    rows.append(dict(source=int(a['source'][i]),path=lookup[int(a['source'][i])],path_over_c_ns=a['x'][i],
        projected_time_residual_ns=res[i],native_fields=values))
result=dict(events=len(norm),norm_above_1ns=int((norm>1).sum()),norm_above_10ns=len(ix),rows=rows,
    interpretation='These remaining four-used-layer tails are retained. Native reconstruction/association checks are still required; do not infer clock wrap or new physics from a stored-time outlier.',
    no_new_selection=True,RICH_read=False,MC_read=False)
r.write(BASE/'used_layer_audit/four_used_tails.json',result)
print(json.dumps(r.native(result),indent=2),flush=True)
