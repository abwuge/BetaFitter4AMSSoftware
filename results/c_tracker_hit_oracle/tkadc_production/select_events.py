import sys
import numpy as np
import uproot
src,dst,ids=sys.argv[1:]
wanted=np.loadtxt(ids,dtype=np.int64,ndmin=2)
with uproot.open(src) as f:
 t=f['amstreea']
 a=t.arrays(['run','event'],library='np')
 keys=np.rec.fromarrays([a['run'].astype(np.int64),a['event'].astype(np.int64)],names='run,event')
 targets=np.rec.fromarrays([wanted[:,0],wanted[:,1]],names='run,event')
 keep=np.isin(keys,targets)
 if keep.sum()!=len(wanted):
  raise RuntimeError('Matched %d of %d target events'%(keep.sum(),len(wanted)))
 arrays={}
 for name in t.keys():
  arrays[name]=t[name].array(library='np')[keep]
 for name in ['tk_adc','tk_qlqs','tk_stripx']:
  if name not in arrays: raise RuntimeError('Missing '+name)
 with uproot.recreate(dst) as out:
  out['amstreea']=arrays
with uproot.open(dst) as f:
 assert f['amstreea'].num_entries==len(wanted)
print('SELECTED',len(wanted),'events',len(arrays),'branches')
