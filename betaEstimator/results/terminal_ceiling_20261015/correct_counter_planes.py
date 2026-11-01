#!/usr/bin/env python3
"""Use actual counter-center planes, not TOFGeom mean-layer planes."""
import re,hashlib
from pathlib import Path
import numpy as np
from audit_common import HERE,C,write_json,metric
SOURCE=Path('/nas01/home/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/CC/Tofdbc.C')
text=SOURCE.read_text();section=text.split('const float TOFGeom::Sci_pz[NLAY][NBAR]=',1)[1].split('};',1)[0]
rows=re.findall(r'\{([^{}]*)\}',section)
coords=np.array([[float(x) for x in row.split(',') if x.strip()] for row in rows],float)
assert coords.shape==(4,10)
b=np.load(str(HERE/'broad_Z6.npz'));old=np.load(str(HERE/'native_collected.npz'));a={k:old[k] for k in old.files}
a['z_counter']=coords[np.arange(4)[None,:],b['tof_barid'].astype(int)]
a['tof_time_plane']=a['tof_time']+(a['z_counter']-a['tof_pos'][:,:,2])/(C*a['tof_dir_z']*a['tof_beta'])
a['primary_four']&=np.all(a['z_counter']!=0,axis=1)&np.isfinite(a['tof_time_plane']).all(1)
a['primary_four_unique']&=a['primary_four'];a['primary_postS4']&=a['primary_four']
x=b['tof_leng'].astype(float)/C;xc=x-x.mean(1)[:,None];w=xc/np.sum(xc*xc,axis=1)[:,None]
a['u_time_plane']=np.sum(w*a['tof_time_plane'],axis=1);a['delta_prop']=a['u_time_plane']-a['u_true']
a['noise_vector']=b['tof_tl'].astype(float)-a['tof_time_plane'];a['noise_vector']-=a['noise_vector'].mean(1)[:,None]
a['noise_slope']=np.sum(w*a['noise_vector'],axis=1);a['noise_orthogonal']=a['noise_vector']-a['noise_slope'][:,None]*xc
valid=a['primary_four'];closure=b['u']-a['u_true']-a['delta_prop']-a['noise_slope']
assert np.max(np.abs(closure[valid]))<1e-12
np.savez_compressed(str(HERE/'native_counter_collected.npz'),**a)
write_json(HERE/'counter_plane_provenance.json',dict(source=str(SOURCE),source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),counter_z=coords,ntuple_path_definition='analysis_amsd69n.C:3424-3429 uses tofclh->Coo as Track interpolate target; tof_pos at lines3387-3390 uses TOFGeom::GetMeanZ and is NOT the path endpoint.',primary_four=int(valid.sum()),delta_prop=metric(a['delta_prop'][valid]),noise=metric(a['noise_slope'][valid]),mean_layer_convention_minus_counter=metric((old['u_time_plane']-a['u_time_plane'])[valid]),max_algebraic_closure=float(np.max(np.abs(closure[valid])))))
print('COUNTER_PLANE_COMPLETE',int(valid.sum()),metric(a['delta_prop'][valid]),metric(a['noise_slope'][valid]),flush=True)
