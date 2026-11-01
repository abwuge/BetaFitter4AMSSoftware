#!/usr/bin/env python3
"""One-way final truth evaluation; refuses changed frozen artifacts."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np,uproot
import reference_model as ref
import identifiable_model as ident
import timing_model as old
from extract_inputs import load_raw,parse_sources,CHARGES
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def predict_all(d,z,frozen,models,legacy):
    prod=1/d['tof_betah'].astype(float)
    ok=old.clean(d);ids=np.flatnonzero(ok)
    result={name:prod.copy() for name in ['ols','old_nll','identifiable']+list(models[z])+['selected']}
    for start in range(0,len(ids),30000):
        ii=ids[start:start+30000];p={k:v[ii] for k,v in d.items()}
        u0=old.kinematics(p)[2];result['ols'][ii]=u0
        result['old_nll'][ii]=old.predict(p,legacy[z])[0]
        result['identifiable'][ii]=ident.predict(p,legacy[z])[0]
        for name,m in models[z].items():result[name][ii]=ref.predict(p,m)[0]
    chosen=frozen['charges'][str(z)]['selected'];name=chosen['name']
    if name=='old_nll':result['selected']=result[name].copy()
    else:result['selected']=result['ols']+chosen['shrink']*(result[name]-result['ols'])
    for name,u in result.items():
        good=np.isfinite(u)&(u>0)
        u[~good]=prod[~good]
    result['producer']=prod
    return result,ok
def describe(e):
    q=np.percentile(e,[16,50,84]);w=(q[2]-q[0])/2
    return dict(n=len(e),median=float(q[1]),width68=float(w),rms=float(np.sqrt(np.mean(e*e))),
                tail_3width=float(np.mean(np.abs(e-q[1])>3*w)))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--models',default='results/reference_training')
    ap.add_argument('--output',default='results/reference_comparison');args=ap.parse_args()
    root=Path(args.models);dest=Path(args.output)
    frozen=json.loads((root/'frozen_summary.json').read_text())
    for name,h in frozen['code_sha256'].items():
        if digest(Path('analysis')/name)!=h:raise ValueError('Changed frozen code: '+name)
    for name,h in frozen['model_sha256'].items():
        if digest(root/'models'/name)!=h:raise ValueError('Changed model: '+name)
    if digest('config/reference_mc_holdout.conf')!=frozen['fresh_sources_sha256']:raise ValueError('Changed holdout')
    models={};legacy={}
    for z in CHARGES:
        models[z]={p.stem.split('_',1)[1]:json.loads(p.read_text()) for p in (root/'models').glob('Z%d_*.json'%z)}
        p=Path('results/data_training/models')/('model_Z%d.json'%z)
        if digest(p)!=frozen['legacy_models_sha256'][str(z)]:raise ValueError('Changed legacy model')
        legacy[z]=json.loads(p.read_text())
    dest.mkdir(parents=True,exist_ok=True)
    audit={'model_frozen_utc':frozen['frozen_utc'],'truth_opened_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
           'training_summary_sha256':digest(root/'frozen_summary.json'),'sources':[],'samples':{}}
    old_ids={z:set() for z in CHARGES}
    for sample,manifest in [('historical','config/sources.conf'),('fresh','config/reference_mc_holdout.conf')]:
        stores={z:[] for z in CHARGES};iv=0
        for sid,(role,z,path) in enumerate(parse_sources(manifest),1):
            if role!='validation':continue
            iv+=1
            d,_=load_raw(path)
            pm={}
            for zz in CHARGES:
                mask=d['charge']==zz
                p,ok=predict_all({k:v[mask] for k,v in d.items()},zz,frozen,models,legacy)
                for name,u in p.items():
                    if name not in pm:pm[name]=np.zeros(len(mask),dtype=float)
                    pm[name][mask]=u
            # Every detector prediction for this source is fixed before truth access.
            with uproot.open(path) as file:
                truth=file['amstreea'].arrays(['mch','mmass','mevmom1','mtrmom','mtrpri','mtrcoo'],library='np')
            truth={k:v[d['source_entry']] for k,v in truth.items()}
            mass=truth['mmass'].astype(float);p=truth['mevmom1'][:,10].astype(float)
            with np.errstate(invalid='ignore',divide='ignore'):
                beta=p/np.sqrt(p*p+mass*mass)
                pl5=truth['mtrmom'][:,4].astype(float)
                l5=pl5/np.sqrt(pl5*pl5+mass*mass)
            l5valid=(truth['mtrpri'][:,4]==1)&(pl5>0)&np.isfinite(l5)&(np.abs(truth['mtrcoo'][:,4,2]-1.698)<.03)
            historical=(truth['mevmom1'][:,10]!=-1000)&(mass>0)
            if sample=='historical':
                part=Path('../betaFitter/results/hybrid_physical_20260913/validation/parts')/('part_%03d_Z%d.root'%(iv,z))
                with uproot.open(str(part)) as file:b=file['betaTree'].arrays(library='np')
                if len(b['Z'])!=int(historical.sum()):raise ValueError('Source size mismatch')
                if not np.array_equal(b['Z'],d['charge'][historical]):raise ValueError('Charge alignment')
                if not np.array_equal(b['mcCharge'],truth['mch'][historical].astype(int)):raise ValueError('Truth charge alignment')
                if not np.allclose(b['linearBeta'],d['tof_betah'][historical],rtol=0,atol=0,equal_nan=True):raise ValueError('Timing alignment')
                if not np.allclose(b['mcBeta'],beta[historical],rtol=1e-6,atol=1e-7):raise ValueError('Truth alignment')
                inds=np.flatnonzero(historical)
                common=(b['Z']==z)&(b['mcCharge']==z)&(b['mcBeta']>=.25)&(b['mcBeta']<=1)
                for k in ('linearBeta','nonlinearBeta','hybridBeta'):common&=np.isfinite(b[k])&(b[k]>0)
                inds=inds[common];baseline=b['hybridBeta'][common]
                truthbeta=b['mcBeta'][common].astype(float)
            else:
                common=(d['charge']==z)&(truth['mch']==z)&l5valid&(mass>0)&(l5>=.25)&(l5<=1)
                common&=np.isfinite(d['tof_betah'])&(d['tof_betah']>0)
                inds=np.flatnonzero(common);baseline=d['tof_betah'][inds].astype(float);truthbeta=l5[inds]
            ids=list(zip(d['run'][inds].tolist(),d['event'][inds].tolist()))
            if sample=='historical':old_ids[z].update(ids)
            elif any(x in old_ids[z] for x in ids):raise ValueError('Fresh/historical event overlap')
            row={'run':d['run'][inds],'event':d['event'][inds],'sourceId':np.full(len(inds),sid,dtype=np.int32),
                 'Z':np.full(len(inds),z,dtype=np.int32),'mcBeta':truthbeta,'mcBetaL5':l5[inds],
                 'l5Valid':l5valid[inds],'linearBeta':d['tof_betah'][inds].astype(float),
                 'hybridBeta':baseline,'newBeta':1/pm['selected'][inds]}
            row.update({name+'U':u[inds] for name,u in pm.items()})
            stores[z].append(row)
            audit['sources'].append({'sample':sample,'path':path,'source_id':sid,'common':len(inds),
             'raw_selected':len(d['charge']),'primary_l5':int(l5valid.sum())})
            if iv%16==0:print('%s final truth %d sources'%(sample,iv),flush=True)
        od=dest/sample;od.mkdir(exist_ok=True)
        audit['samples'][sample]={}
        for z,parts in stores.items():
            rows={k:np.concatenate([p[k] for p in parts]) for k in parts[0]}
            ids=np.rec.fromarrays([rows['run'],rows['event']])
            if len(np.unique(ids))!=len(ids):raise ValueError('Duplicate final events')
            with uproot.recreate(str(od/('comparison_Z%d.root'%z))) as f:f['comparison']=rows
            metrics={}
            for low,high in [(.6,.9),(.9,1.000001),(.6,1.000001),(.95,1.000001)]:
                m=(rows['mcBeta']>=low)&(rows['mcBeta']<high);trueu=1/rows['mcBeta'][m]
                if m.sum()<20:continue
                group={'beta_low':low,'beta_high':min(high,1.),'N':int(m.sum()),'methods':{}}
                group['methods']['baseline']=describe(1/rows['hybridBeta'][m]-trueu)
                for name in pm:group['methods'][name]=describe(rows[name+'U'][m]-trueu)
                metrics[str(low)+'_'+str(high)]=group
            audit['samples'][sample][str(z)]={'N':len(rows['run']),'metrics':metrics}
            print('%s Z%d common=%d'%(sample,z,len(rows['run'])),flush=True)
        (od/'evaluation_summary.json').write_text(json.dumps(audit['samples'][sample],indent=2))
    (dest/'evaluation_audit.json').write_text(json.dumps(audit,indent=2))
    print('FINAL EVALUATION COMPLETE: frozen models unchanged',flush=True)
if __name__=='__main__':main()
