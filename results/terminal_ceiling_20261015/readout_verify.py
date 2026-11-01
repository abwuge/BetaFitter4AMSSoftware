import numpy as np,uproot
from audit_common import HERE,C,metric,write_json
b=np.load(str(HERE/'broad_Z6.npz'));n=np.load(str(HERE/'native_counter_collected.npz'));r=uproot.open(str(HERE/'native_readout_048.root'))['readout'].arrays(library='np')
keys=b['run'].astype(np.uint64)*np.uint64(2**32)+b['event'].astype(np.uint64);order=np.argsort(keys);ids=r['run'].astype(np.uint64)*np.uint64(2**32)+r['event'].astype(np.uint64);ix=order[np.searchsorted(keys[order],ids)];assert np.array_equal(keys[ix],ids)
x=b['tof_leng'][ix].astype(float)/C;xc=x-x.mean(1)[:,None];w=xc/np.sum(xc*xc,axis=1)[:,None]
shift=b['tof_tl'][ix]-r['native_time'];shift-=shift.mean(1)[:,None];good=n['truth_valid'][ix]
ans=dict(events=len(ix),found=int(r['found'].sum()),associated=int(r['associated'].sum()),position_finite=int(np.isfinite(r['position']).all((1,2)).sum()),max_counter_z_difference=float(np.nanmax(np.abs(r['position'][:,:,2]-r['counter_z']))),max_counter_table_difference=float(np.nanmax(np.abs(r['counter_z']-n['z_counter'][ix]))),ntuple_minus_native_relative_time_ns={str(j):metric(shift[:,j]) for j in range(4)},ntuple_minus_native_slope=metric(np.sum(w*shift,axis=1)),native_time_error=metric((np.sum(w*r['native_time'],axis=1)-n['u_true'][ix])[good]),ntuple_time_error=metric((b['u'][ix]-n['u_true'][ix])[good]))
if 'tof_rtime' in b:ans['max_raw_side_time_difference']=float(np.nanmax(np.abs(b['tof_rtime'][ix]-r['raw_time'])))
write_json(HERE/'native_readout_verification.json',ans);print(ans,flush=True)
