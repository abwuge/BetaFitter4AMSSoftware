#!/usr/bin/env python3
"""Select and freeze TOF corrections using independent flight measurements only."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
import uproot
import reference_model as ref
import identifiable_model as ident
import timing_model as old
from extract_inputs import CHARGES
from extract_reference import ALL_BRANCHES
from augment_reference import FIELDS
BINS=np.array([.78,.85,.90,.95,.97,.985,1.02])
ALPHAS=(.01,.1,1.)
SHRINK=(0.,.25,.5,.75,1.)
CACHE_KEYS=set(['run','event','tof_betah','tof_tl','tof_leng','tof_barid','tof_edep','tof_bpatt',
 'tk_edep','position','direction','charge','source_entry','source_id'])|set(ALL_BRANCHES[14:])|set(FIELDS)
def subset(d,mask):return {k:v[mask] for k,v in d.items()}
def load(path):
    with np.load(str(path),allow_pickle=False) as f:d={k:f[k] for k in f.files}
    expected=set(['run','event','tof_betah','tof_tl','tof_leng','tof_barid','tof_edep','tof_bpatt','tk_edep','position','direction','charge','source_entry','source_id']) | set(__import__('extract_reference').REFERENCE_BRANCHES) | set(FIELDS)
    forbidden=sorted(set(d)-expected)
    if set(d)!=expected:raise ValueError('Cache schema differs from measured allowlist: '+str(sorted(set(d)^expected)))
    if forbidden:raise ValueError('Truth-bearing reference cache: '+str(forbidden))
    return d
def quality(d):
    b=d['rich_beta'][:,0].astype(float);naf=d['rich_NaF'].astype(bool)
    pos=d['rich_pos'];xy=np.max(np.abs(pos[:,:2]),axis=1);rad=np.sum(pos[:,:2]**2,axis=1)
    npe=d['rich_npe'];bc=d['rich_BetaConsistency'];pc=d['rich_PMTChargeConsistency']
    ok=(d['irich']>=0)&(d['itrtrack']>=0)&(d['rich_itrtrack']==d['itrtrack'])
    ok&=d['rich_good']&d['rich_clean']&(d['rich_pb']>=.01)&(d['rich_pmt']>=3)&(d['rich_used']>=4)
    ok&=np.isfinite(b)&(b>=.78)&(b<=1.02)&np.where(naf,b>=.78,b>=.955)
    ok&=np.isfinite(d['rich_ebeta1'])&(d['rich_ebeta1']>0)&(d['rich_ebeta1']<.005)
    ok&=np.isfinite(bc)&(bc>=0)&(bc<=np.where(naf,.01,.005))
    ok&=np.isfinite(pc)&(pc>=0)&(pc<=15)&np.all(np.isfinite(npe),axis=1)&(npe[:,2]>0)
    ok&=(npe[:,0]>=.4*npe[:,2])&(npe[:,1]>=np.where(naf,1.,2.))
    geo=np.where(naf,xy<=17,(xy>=19)&(rad<=3500)&~np.isin(d['rich_tile'],[3,7,12,20,87,100,108]))
    ok&=geo&np.all(np.isfinite(pos),axis=1)&old.clean(d)
    return ok
def inner_fold(runs):
    table={int(r):int(hashlib.sha256(('reference-inner:'+str(int(r))).encode()).hexdigest()[:8],16)%4==0 for r in np.unique(runs)}
    return np.array([table[int(r)] for r in runs])
def summary(u,b,base=None):
    e=u-1/b;out=[]
    for lo,hi in zip(BINS[:-1],BINS[1:]):
        m=(b>=lo)&(b<hi if hi<1.02 else b<=hi)
        if m.sum()<40:continue
        v=e[m];q=np.percentile(v,[16,50,84]);width=float((q[2]-q[0])/2)
        row={'beta_low':float(lo),'beta_high':float(hi),'n':int(m.sum()),
             'median':float(q[1]),'width68':width,'rms':float(np.sqrt(np.mean(v*v)))}
        if base is not None:
            bv=(base-1/b)[m];bq=np.percentile(bv,[16,50,84]);bw=(bq[2]-bq[0])/2
            row['width_ratio']=width/max(bw,1e-10)
            # Same baseline scale in all candidates; reference risk, no truth.
            delta=3*max(bw,1e-5)
            def huber(x):
                a=np.abs(x);return np.mean(np.where(a<=delta,.5*x*x,delta*(a-.5*delta)))
            row['huber_risk_ratio']=float(huber(v)/max(huber(bv),1e-16))
        out.append(row)
    return out
def score(rows):
    usable=[r for r in rows if r['n']>=200]
    if not usable:return float('inf')
    # Direct velocity residual objective, balancing reference-beta intervals.
    ratios=np.array([.5*r['width_ratio']**2+.5*r['huber_risk_ratio'] for r in usable])
    return float(.75*np.mean(ratios)+.25*np.max(ratios))
def capped(d,limit,seed):
    if len(d['run'])<=limit:return d
    idx=np.sort(np.random.RandomState(seed).choice(len(d['run']),limit,replace=False))
    return subset(d,idx)
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cache',default='/tmp/betaEstimator-reference')
    ap.add_argument('--output',default='results/reference_training');ap.add_argument('--max-training',type=int,default=100000)
    args=ap.parse_args();dest=Path(args.output);(dest/'models').mkdir(parents=True,exist_ok=True)
    audit=json.loads((Path(args.cache)/'reference_audit.json').read_text())
    report={'truth_access':False,'reference':'independent measured RICH; high-rigidity anchor cross-check',
      'audit':audit,'charges':{},'alphas':ALPHAS,'shrinkage':SHRINK,'selection':'flight inner-run validation only',
      'reference_plane':'RICH radiator; calibration restricted to beta_RICH>=0.99 to suppress physical center-to-radiator speed change; no exact transport assumed'}
    for z in CHARGES:
        d=load(Path(args.cache)/('train_Z%d.npz'%z));v=load(Path(args.cache)/('validation_Z%d.npz'%z))
        if set(d['run'])&set(v['run']):raise ValueError('Outer run overlap')
        inn=inner_fold(d['run']);q=quality(d);vq=quality(v)
        train=subset(d,q&~inn&(d['rich_beta'][:,0]>=.99));cv=subset(d,q&inn&(d['rich_beta'][:,0]>=.99));va=subset(v,vq)
        if len(train['run'])<300 or len(cv['run'])<100:raise ValueError('Insufficient independent reference for Z%d'%z)
        train=capped(train,args.max_training,813+z)
        _,_,u0=old.kinematics(train);target=u0-1/train['rich_beta'][:,0]
        _,_,cb=old.kinematics(cv);rb=cv['rich_beta'][:,0].astype(float)
        legacy=json.loads((Path('results/data_training/models')/('model_Z%d.json'%z)).read_text())
        records=[];candidates={}
        for mode in ('geometry','full'):
            for alpha in ALPHAS:
                name='rich_%s_a%s'%(mode,alpha)
                model=ref.fit(train,z,target,alpha,mode=mode)
                _,corr=ref.predict(cv,model)
                candidates[name]=model
                for shrink in SHRINK:
                    rows=summary(cb-shrink*corr,rb,cb)
                    records.append({'name':name,'shrink':shrink,'score':score(rows),'bins':rows})
                print('Z%d fitted %s N=%d'%(z,name,len(train['run'])),flush=True)
        # High-rigidity calibration is not an inference input or a mass hypothesis.
        anchor=subset(d,(~inn)&np.isfinite(d['tk_rigidity'])&(d['tk_rigidity']>100)&old.clean(d))
        if len(anchor['run'])>=300:
            anchor=capped(anchor,args.max_training,1813+z)
            _,_,au=old.kinematics(anchor)
            am=ref.fit(anchor,z,au-1.,.1,mode='full');name='saturation_full'
            candidates[name]=am;_,corr=ref.predict(cv,am)
            for shrink in SHRINK:
                rows=summary(cb-shrink*corr,rb,cb)
                records.append({'name':name,'shrink':shrink,'score':score(rows),'bins':rows})
        iu,ic=ident.predict(cv,legacy)
        for shrink in SHRINK:
            rows=summary(cb-shrink*ic,rb,cb)
            records.append({'name':'identifiable','shrink':shrink,'score':score(rows),'bins':rows})
        ou=old.predict(cv,legacy)[0]
        rows=summary(ou,rb,cb)
        records.append({'name':'old_nll','shrink':1.,'score':score(rows),'bins':rows})
        # Tie breaks prefer stronger ridge and less applied correction.
        best=min(records,key=lambda r:(r['score'],r['shrink'],r['name']))
        for name,model in candidates.items():
            p=dest/'models'/('Z%d_%s.json'%(z,name));p.write_text(json.dumps(model,indent=2,allow_nan=False))
        # Models trained on inner-training runs stay fixed; no outer-holdout refit.
        _,_,vb=old.kinematics(va);vr=va['rich_beta'][:,0].astype(float)
        preds={'linear':vb,'old_nll':old.predict(va,legacy)[0],'identifiable':ident.predict(va,legacy)[0]}
        for name,model in candidates.items():preds[name]=ref.predict(va,model)[0]
        if best['name']=='old_nll':pred=preds['old_nll']
        else:pred=vb+best['shrink']*(preds[best['name']]-vb)
        metrics={name:summary(u,vr,vb) for name,u in preds.items()}
        metrics['selected']=summary(pred,vr,vb)
        with uproot.recreate(str(dest/('flight_Z%d.root'%z))) as out:
            rows={'run':va['run'],'event':va['event'],'sourceId':va['source_id'],'richBeta':vr,
                  'linearU':vb,'oldNllU':preds['old_nll'],'identifiedU':preds['identifiable'],'selectedU':pred,
                  'NaF':va['rich_NaF']}
            out['reference']=rows
        report['charges'][str(z)]={'train_reference':len(train['run']),'inner_reference':len(cv['run']),
            'holdout_reference':len(va['run']),'anchor_train':len(anchor['run']),
            'train_runs':len(np.unique(train['run'])),'inner_runs':len(np.unique(cv['run'])),
            'holdout_runs':len(np.unique(va['run'])),'selected':best,'inner_candidates':records,
            'holdout_metrics':metrics,'coverage':{'train_NaF':int(train['rich_NaF'].sum()),
               'holdout_NaF':int(va['rich_NaF'].sum())}}
        print('Z%d selected %s shrink %.2f inner %.5f holdout N=%d'%(z,best['name'],best['shrink'],best['score'],len(va['run'])),flush=True)
    report['frozen_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
    names=['reference_model.py','identifiable_model.py','timing_model.py','extract_reference.py','augment_reference.py','train_reference.py','extract_inputs.py','evaluate_reference.py']
    report['code_sha256']={n:digest(Path('analysis')/n) for n in names}
    report['model_sha256']={p.name:digest(p) for p in sorted((dest/'models').glob('*.json'))}
    report['legacy_models_sha256']={str(z):digest(Path('results/data_training/models')/('model_Z%d.json'%z)) for z in CHARGES}
    report['fresh_sources_sha256']=digest('config/reference_mc_holdout.conf')
    (dest/'frozen_summary.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    print('FROZEN: all models and selections; no MC truth read.',flush=True)
if __name__=='__main__':main()
