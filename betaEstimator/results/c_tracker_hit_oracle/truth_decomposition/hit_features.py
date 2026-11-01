#!/usr/bin/env python3
import numpy as np,uproot,json
from pathlib import Path
H=Path(__file__).resolve().parent;a=np.load(str(H/'aligned.npz'))
ids=[7,8,9,10,11,12,13];c=a['mevcoo1'][:,ids,:].astype(float);d=a['mevdir1'][:,ids,:].astype(float);h=a['tk_hitc'][:,1:8,:].astype(float)
t=c+d/d[:,:,2,None]*(h[:,:,2]-c[:,:,2])[:,:,None]
y=(h[:,:,1]-t[:,:,1])*1e4
ok=(a['mtrpri'][:,1:8]==1)&(a['mevmom1'][:,ids]>0)&np.isfinite(y)
adc=a['tk_adc'][:,1:8,:,:].astype(float);adc=np.maximum(adc,0);total=adc.sum(axis=-1);frac=adc/np.maximum(total[:,:,:,None],1e-6)
role=np.where(a['ref_role']==b'train',0,np.where(a['ref_role']==b'valid',1,2)).astype('int32')
N=len(role);columns=[];names=[]
def add(v,name):
 v=np.asarray(v).reshape(N,7,-1)
 for j in range(v.shape[-1]):columns.append(v[:,:,j].reshape(-1));names.append(name+'_'+str(j))
add(np.broadcast_to(np.arange(2,9),(N,7)),'layer');add(a['tk_lid'][:,1:8],'tkid');add(a['tk_cof'][:,1:8,:],'cog');add(a['tk_dir'][:,1:8,:2]/a['tk_dir'][:,1:8,2,None],'slope');add(h[:,:,:2],'globalxy')
base=len(names)
add(frac,'adc_fraction');add(np.log(np.maximum(total,1e-6)),'log_adc_total')
X=np.column_stack(columns).astype('float32');keep=ok.ravel()&np.isfinite(X).all(axis=1)
payload={'hit_index':np.arange(N*7,dtype='int32')[keep],'role':np.repeat(role,7)[keep],'event_index':np.repeat(a['ref_event_index'],7).astype('int32')[keep],'source':np.repeat(a['ref_source'],7).astype('int32')[keep],'target':np.clip(y.ravel()[keep],-60,60).astype('float32'),'true_error_um':y.ravel()[keep].astype('float32')}
for j in range(X.shape[1]):payload['f%02d'%j]=X[keep,j]
with uproot.recreate(str(H/'hit_features.root')) as f:
 f.mktree('events',{k:v.dtype for k,v in payload.items()});f['events'].extend(payload)
(H/'hit_features.json').write_text(json.dumps({'features':names,'base_count':base,'hits':int(sum(keep)),'definition':'local measured hit error target; absolute ADC amplitude allowed only to infer hit response, no R/beta/isotope inputs'},indent=2))
for label,count in [('base',base),('adc',len(names))]:
 (H/('hit_'+label+'.txt')).write_text('\n'.join('f%02d'%j for j in range(count))+'\n')
print('HITS',sum(keep),'features',len(names),'base',base)
