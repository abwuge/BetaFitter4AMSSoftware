#!/usr/bin/env python3
import numpy as np, uproot, json
from pathlib import Path
H=Path(__file__).resolve().parent
A=np.load(str(H/'aligned.npz'))
with uproot.open(str(H/'native.root')) as f:
    N=f['probe'].arrays(['event_index','idealrec','jacrec'],library='np')
with uproot.open(str(H/'reverse.root')) as f:
    E=f['reverse'].arrays(['event_index','yref','jac'],library='np')
assert np.array_equal(N['event_index'],A['ref_event_index'])
assert np.array_equal(E['event_index'],A['ref_event_index'])
y=A['tk_hitc'][:,1:8,1].astype(np.float64)
z=A['tk_hitc'][:,1:8,2].astype(np.float64)
refN=N['idealrec'][:,1:8,1].astype(np.float64)
jacN=N['jacrec'][:,1:8,1,:][:,:,[1,3,4]].astype(np.float64)
refE=E['yref'].astype(np.float64)
jacE=E['jac'].astype(np.float64)
roles=A['ref_role']
train=roles==b'train'

def batch_fit(J,r,idx,npar):
    jj=J[:,idx,:npar]
    rr=r[:,idx]
    ata=np.einsum('nki,nkj->nij',jj,jj)
    # tiny relative ridge only for numerical stability
    tr=np.trace(ata,axis1=1,axis2=2)
    ata += np.eye(npar)[None,:,:]*(1e-12+1e-12*tr[:,None,None])
    atr=np.einsum('nki,nk->ni',jj,rr)
    return np.linalg.solve(ata,atr)

def engineer(ref,J):
    r=y-ref
    dz=np.diff(z,axis=1)
    slopes=np.diff(r,axis=1)/dz
    kinks=np.diff(slopes,axis=1)
    cols=[r,slopes,kinks]
    names=(['res_L%d'%l for l in range(2,9)] +
           ['resSlope_L%dL%d'%(l,l+1) for l in range(2,8)] +
           ['resKink_L%d'%l for l in range(3,8)])
    # Global secondary fit: how much the measured points want to move the R-conditioned reference.
    c2=batch_fit(J,r,list(range(7)),2)
    c3=batch_fit(J,r,list(range(7)),3)
    cols += [c2,c3]
    names += ['global_holdR_dy','global_holdR_dty','global_freeR_dy','global_freeR_dty','global_freeR_dlogR']
    # Leave-one-layer cross predictions.
    loo_hold=[]; loo_free=[]; loo_dR=[]
    for i in range(7):
        keep=[j for j in range(7) if j!=i]
        b2=batch_fit(J,r,keep,2); b3=batch_fit(J,r,keep,3)
        p2=np.einsum('ni,ni->n',J[:,i,:2],b2)
        p3=np.einsum('ni,ni->n',J[:,i,:3],b3)
        loo_hold.append(r[:,i]-p2); loo_free.append(r[:,i]-p3); loo_dR.append(b3[:,2])
    cols += [np.column_stack(loo_hold),np.column_stack(loo_free),np.column_stack(loo_dR)]
    names += ['looHold_L%d'%l for l in range(2,9)]
    names += ['looFree_L%d'%l for l in range(2,9)]
    names += ['looOther_dlogR_L%d'%l for l in range(2,9)]
    # Leave one adjacent pair out; this is a stronger anti-circularity test for a local region.
    pair_hold=[]; pair_free=[]; pair_dR=[]; pair_norm=[]
    for i in range(6):
        omit=[i,i+1]; keep=[j for j in range(7) if j not in omit]
        b2=batch_fit(J,r,keep,2); b3=batch_fit(J,r,keep,3)
        qh=[]; qf=[]
        for j in omit:
            qh.append(r[:,j]-np.einsum('ni,ni->n',J[:,j,:2],b2))
            qf.append(r[:,j]-np.einsum('ni,ni->n',J[:,j,:3],b3))
        qh=np.column_stack(qh); qf=np.column_stack(qf)
        pair_hold += [qh[:,0],qh[:,1]]
        pair_free += [qf[:,0],qf[:,1]]
        pair_dR.append(b3[:,2])
        pair_norm.append(np.sqrt(np.sum(qf*qf,axis=1)))
    cols += [np.column_stack(pair_hold),np.column_stack(pair_free),np.column_stack(pair_dR),np.column_stack(pair_norm)]
    for kind in ['pairHold','pairFree']:
        names += ['%s_L%dL%d_end%d'%(kind,l,l+1,e) for l in range(2,8) for e in (0,1)]
    names += ['pairOther_dlogR_L%dL%d'%(l,l+1) for l in range(2,8)]
    names += ['pairFreeNorm_L%dL%d'%(l,l+1) for l in range(2,8)]
    # Independent upper/lower curvature preferences.
    top=batch_fit(J,r,[0,1,2,3],3)
    bot=batch_fit(J,r,[3,4,5,6],3)
    cols += [top,bot,(top-bot)]
    names += ['top_dy','top_dty','top_dlogR','bottom_dy','bottom_dty','bottom_dlogR',
              'topMinusBottom_dy','topMinusBottom_dty','topMinusBottom_dlogR']
    # Event summaries emphasize visible-kink / large-scatter cases.
    summary=np.column_stack([
        np.sqrt(np.mean(r*r,axis=1)),np.max(np.abs(r),axis=1),
        np.sqrt(np.mean(slopes*slopes,axis=1)),np.max(np.abs(slopes),axis=1),
        np.sqrt(np.mean(kinks*kinks,axis=1)),np.max(np.abs(kinks),axis=1)
    ])
    cols.append(summary)
    names += ['res_rms','res_absmax','slope_rms','slope_absmax','kink_rms','kink_absmax']
    X=np.column_stack(cols)
    assert X.shape[1]==len(names),(X.shape,len(names))
    return X,names

XN,namesN=engineer(refN,jacN)
XE,namesE=engineer(refE,jacE)
assert namesN==namesE
# Train-only robust preprocessing, separately for N/E so neither test nor validation sets scale.
def prep(X):
    lo=np.nanpercentile(X[train],0.5,axis=0); hi=np.nanpercentile(X[train],99.5,axis=0)
    T=np.minimum(np.maximum(X[train],lo),hi)
    med=np.nanmedian(T,axis=0); q25,q75=np.nanpercentile(T,[25,75],axis=0)
    scale=np.maximum(0.741301109*(q75-q25),0.1*np.nanstd(T,axis=0))
    good=np.isfinite(med)&np.isfinite(scale)&(scale>1e-12)
    Y=np.minimum(np.maximum(X[:,good],lo[good]),hi[good])
    Y=np.where(np.isfinite(Y),Y,med[good])
    Y=np.clip((Y-med[good])/scale[good],-8,8).astype('float32')
    return Y,[namesN[i] for i in np.flatnonzero(good)]
ZN,nN=prep(XN); ZE,nE=prep(XE)
# Align existing safe baseline feature tree by global event_index.
base_path=H.parent/'tkadc_shape_features.root'
base_list=[x.strip() for x in (H.parent/'features_tkadc_base.txt').read_text().splitlines() if x.strip()]
with uproot.open(str(base_path)) as f:
    B=f['events'].arrays(['event_index','role','source','bin','R','R_true','u','u_true','target']+base_list,library='np')
lookup={int(v):i for i,v in enumerate(A['ref_event_index'])}
ix=np.asarray([lookup[int(v)] for v in B['event_index']],dtype=np.int64)
assert np.array_equal(B['role'],np.where(A['ref_role'][ix]==b'train',0,np.where(A['ref_role'][ix]==b'valid',1,2)))
payload={k:B[k] for k in ['event_index','role','source','bin','R','R_true','u','u_true','target']}
payload['label']=np.zeros(len(ix),dtype=np.int8)
for b in base_list: payload[b]=B[b]
for j in range(ZN.shape[1]): payload['n%03d'%j]=ZN[ix,j]
for j in range(ZE.shape[1]): payload['e%03d'%j]=ZE[ix,j]
out=H/'reverse_features.root'
with uproot.recreate(str(out),compression=uproot.ZLIB(4)) as f: f['events']=payload
(H/'features_reverse_noeloss.txt').write_text('\n'.join('n%03d'%j for j in range(ZN.shape[1]))+'\n')
(H/'features_reverse_eloss.txt').write_text('\n'.join('e%03d'%j for j in range(ZE.shape[1]))+'\n')
(H/'features_reverse_base.txt').write_text('\n'.join(base_list)+'\n')
(H/'features_reverse_base_eloss.txt').write_text('\n'.join(base_list+['e%03d'%j for j in range(ZE.shape[1])])+'\n')
# Energy-loss trajectory impact, measured in microns.
diff=(refE-refN)*1e4
q=np.percentile(diff,[16,50,84],axis=0)
meta={'feature_count_noeloss':int(ZN.shape[1]),'feature_count_eloss':int(ZE.shape[1]),
      'base_count':len(base_list),'names':nE,
      'mean_eloss_reference_shift_um':[
        {'layer':l,'median':float(q[1,i]),'w68':float((q[2,i]-q[0,i])/2)}
        for i,l in enumerate(range(2,9))]}
(H/'reverse_features.json').write_text(json.dumps(meta,indent=2))
print(json.dumps(meta,indent=2))
