import numpy as np
from pathlib import Path
a=np.load(str(Path(__file__).with_name('aligned.npz')))
m=a['ref_role']==b'test'; ids=np.array([0,7,8,9,10,11,12,13,18])
c=a['mevcoo1'][:,ids,:].astype(float); d=a['mevdir1'][:,ids,:].astype(float); h=a['tk_hitc'].astype(float)
t=c+d/d[:,:,2,None]*(h[:,:,2]-c[:,:,2])[:,:,None]
def w(x):return float(np.diff(np.percentile(x,[16,84]))[0]/2)
print('baseline',w(np.log(a['ref_R'][m]/a['ref_R_true'][m])))
for l in range(9):
 good=m & ((a['tk_hitb'][:,1]&(1<<l))!=0)&(c[:,l,2]>-999)&(a['mtrpri'][:,l]==1)
 print('layer',l+1,'N',sum(good),'dz',np.percentile(h[good,l,2]-c[good,l,2],[16,50,84]),'Y h-true projected um',np.percentile((h[good,l,1]-t[good,l,1])*1e4,[1,16,50,84,99]),'Y h-mccluster um',np.percentile((h[good,l,1]-a['mtrcoo1'][good,l,1])*1e4,[16,50,84]),'pdelta',np.percentile(np.log(a['mevmom1'][good,ids[l]]/(6*a['ref_R_true'][good])),[16,50,84]))
print('sample coords',c[m][0],h[m][0]);print('bits',np.unique(a['tk_hitb'][m],axis=0,return_counts=True));print('ptrue exact',np.max(abs(a['mevmom1'][:,10]/6-a['ref_R_true'])))
