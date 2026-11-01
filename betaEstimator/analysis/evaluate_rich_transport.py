#!/usr/bin/env python3
"""Evaluate the frozen center/RICH model on existing MC and flight sources."""
import argparse
import csv
import hashlib
import json
import multiprocessing as mp
import time
from pathlib import Path
import numpy as np
import uproot
import rich_transport as rt

ROOT = rt.ROOT
EXP = ROOT / 'results/c_joint_beta_rigidity_closure/expanded_validation'
MEASURED = ['run','event','tof_tl','tof_leng','tof_pos','tof_edep','tof_barid','tk_ql','tk_rigidityk','tk_pos','tk_dir',
            'irich','itrtrack','rich_itrtrack','rich_beta','rich_good','rich_clean','rich_pb','rich_pmt','rich_used',
            'rich_ebeta1','rich_BetaConsistency','rich_PMTChargeConsistency','rich_npe','rich_NaF','rich_pos','rich_tile']
SELECTION = ['tk_qin','tk_q','tof_bpatt','tk_erigidity1','tk_chis1']
EVAL = ['mevmom1','mevcoo1','mtrmom','mtrpri','mtrcoo']
MC_IDS = None


def digest(path):
    h = hashlib.sha256()
    with open(str(path),'rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()


def dump(path, data):
    def conv(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(data,indent=2,default=conv,allow_nan=True)+'\n')


def sub(d, rows):
    return {k:v[rows] for k,v in d.items()}


def concatenate(pieces):
    return {k:np.concatenate([p[k] for p in pieces]) for k in pieces[0]}


def match_rows(run,event,target_run,target_event):
    ids = np.rec.fromarrays([run.astype(np.int64),event.astype(np.int64)],names='run,event')
    target = np.rec.fromarrays([target_run.astype(np.int64),target_event.astype(np.int64)],names='run,event')
    common,i,j = np.intersect1d(ids,target,return_indices=True)
    if len(common)!=len(target):raise RuntimeError('MC source/event join is incomplete')
    return i[np.argsort(j)]


def mc_file(task):
    label, source, path = task
    rows = np.flatnonzero((MC_IDS['label']==label)&(MC_IDS['source']==source))
    with uproot.open(path) as f:
        t=f['amstreea']
        a=t.arrays(MEASURED,library='np')
        take=match_rows(a['run'],a['event'],MC_IDS['run'][rows],MC_IDS['event'][rows])
        out=sub(a,take)
        # All model files and the protocol have already been frozen before
        # this evaluation-only branch read. These keys never enter prediction.
        e=t.arrays(EVAL,library='np')
        for k in EVAL:out['eval_'+k]=e[k][take]
        raw=int(t.num_entries)
    out['row']=rows
    return out,{'label':label,'source':source,'path':path,'raw':raw,'selected':len(rows)}


def flight_file(task):
    source, role, path = task
    with uproot.open(path) as f:
        t=f['amstreea']
        a=t.arrays(MEASURED+SELECTION,library='np')
        rows=np.flatnonzero(rt.core_selection(a))
        out={k:a[k][rows] for k in MEASURED}
        out['source']=np.full(len(rows),source,np.int32)
        out['outer_holdout']=np.full(len(rows),role=='validation',bool)
        raw=int(t.num_entries)
    return out,{'source':source,'role':role,'path':path,'raw':raw,'selected':len(rows)}


def extract(kind,workers,dest):
    global MC_IDS
    pieces,audit=[],[]
    if kind=='mc':
        with np.load(str(EXP/'predictions.npz')) as f: MC_IDS={k:f[k] for k in f.files}
        tasks=[]
        with (EXP/'manifest.tsv').open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                tasks.append(({'C12':0,'C13':1}[row['isotope']],int(row['source']),row['path']))
        worker=mc_file
    else:
        tasks=[]
        for line in (ROOT/'config/data_sources.conf').read_text().splitlines():
            if line.strip() and not line.startswith('#'):
                role,charge,path=line.split()
                tasks.append((len(tasks),role,path))
        worker=flight_file
    with mp.Pool(workers) as pool:
        for j,(part,record) in enumerate(pool.imap(worker,tasks),1):
            pieces.append(part);audit.append(record)
            if j%24==0 or j==len(tasks):print(kind,'files',j,'/',len(tasks),'selected',sum(r['selected'] for r in audit),flush=True)
    a=concatenate(pieces)
    if kind=='mc':
        order=np.argsort(a['row']);a=sub(a,order)
        if not np.array_equal(a.pop('row'),np.arange(len(MC_IDS['run']))):raise RuntimeError('MC order differs')
        for k in ['run','event']:
            if not np.array_equal(a[k],MC_IDS[k]):raise RuntimeError('MC identity differs')
        for k in ['label','source','u_true','truth_valid']:a['eval_'+k]=MC_IDS[k]
    else:
        ids=np.rec.fromarrays([a['run'],a['event']],names='run,event')
        if len(np.unique(ids))!=len(ids):raise RuntimeError('Duplicate flight events')
    dump(dest/(kind+'_source_audit.json'),audit)
    return a


def width(x):
    return float(np.diff(np.percentile(x,[16,84]))[0]/2)


def stats(x):
    if len(x)==0:return {'n':0}
    return {'n':int(len(x)),'mean':float(np.mean(x)),'median':float(np.median(x)),
            'w68':width(x),'std':float(np.std(x)),'rmse':float(np.sqrt(np.mean(x*x))),
            'p95_abs':float(np.percentile(np.abs(x),95))}


def compare(old,new,reference):
    return {'baseline_beta':stats(old-reference),'new_beta':stats(new-reference),
            'baseline_u':stats(1/old-1/reference),'new_u':stats(1/new-1/reference),
            'gain_beta_pct':100*(1-width(new-reference)/width(old-reference)),
            'gain_u_pct':100*(1-width(1/new-1/reference)/width(1/old-1/reference))}


def bootstrap(old,new,reference,source,label=None,repeats=400):
    if len(old)<20:return {'n':len(old),'repeats':0}
    label=np.zeros(len(old),int) if label is None else label
    groups={}
    for lab in np.unique(label):
        groups[int(lab)]=[np.flatnonzero((label==lab)&(source==s)) for s in np.unique(source[label==lab])]
    rng=np.random.RandomState(20261007)
    eb0=old-reference;eb1=new-reference
    eu0=1/old-1/reference;eu1=1/new-1/reference
    values=[]
    for _ in range(repeats):
        ix=np.concatenate([g[i] for g in groups.values() for i in rng.randint(0,len(g),len(g))])
        values.append([100*(1-width(eb1[ix])/width(eb0[ix])),100*(1-width(eu1[ix])/width(eu0[ix]))])
    values=np.array(values)
    return {'repeats':repeats,'source_groups':{str(k):len(v) for k,v in groups.items()},
            'beta_gain_pct_q025_q50_q975':np.percentile(values[:,0],[2.5,50,97.5]),
            'u_gain_pct_q025_q50_q975':np.percentile(values[:,1],[2.5,50,97.5]),
            'fixed_models':True}


def freeze(dest):
    dest.mkdir(parents=True,exist_ok=True)
    paths=[rt.MODEL,Path(rt.__file__).resolve(),Path(__file__).resolve(),rt.BASE/'run_closure.py',rt.BASE/'run_energy_residualized.py',EXP/'manifest.tsv',ROOT/'config/data_sources.conf']
    hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
    f=dest/'freeze.json'
    if f.exists():
        old=json.loads(f.read_text())
        if old['sha256']!=hashes:raise RuntimeError('Frozen source or model changed; use a separate explicitly labeled study directory')
        return
    info={'utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'sha256':hashes,
          'center_model':'existing 7200-event weak isotope-proportion anchored-energy fit; no retraining',
          'domain':'measured charge 6, 2.5 <= tk_rigidityk[1,1] < 6 GV',
          'transport':'native elemental material map; uncalibrated Bethe mean; midpoint energy integration; cubic Hermite measured tracker path',
          'mass':'continuous reconstructed mass from p=6R and center beta; no isotope label',
          'map_path':rt.MAP,'map_sha256':digest(rt.MAP),
          'RICH_target_z_cm':{'NaF':-75.552,'aerogel':-74.705},
          'radiator_assignment':'tracker extrapolated max(abs(x),abs(y)) at z=-75.1 below 17.5 cm',
          'rich_reference':'rich_beta[0]; quality checks independent of TOF-beta agreement',
          'metrics':'paired W68 beta and inverse-beta residuals; bias/RMSE/tails; source-cluster bootstrap 400 repeats',
          'sensitivity':'material 0.8/1.2 and common mass C12/C13; no selection of best setting',
          'truth_policy':'MC truth read only after this freeze for evaluation; no direct truth checkpoint at RICH',
          'plane_names':rt.PLANE_NAMES,'interval_names':rt.LOSS_NAMES}
    dump(f,info)
    (dest/'center_model.json').write_bytes(rt.MODEL.read_bytes())
    print('FROZEN',f,flush=True)


def evaluate(kind,a,dest):
    model=rt.CenterModel(dest/'center_model.json')
    d=rt.center_features(a)
    unew=model.predict(d)
    out={'run':a['run'],'event':a['event'],'R':d['R'],'bin':d['bin'],'u_raw':d['u'],'u_new':unew,
         'rich_beta':a['rich_beta'][:,0].astype(float),'NaF':a['rich_NaF'].astype(bool),
         'rich_quality':rt.rich_quality(a),'rich_strict':rt.rich_quality(a,True)}
    report={'kind':kind,'selected':len(unew)}
    if kind=='mc':
        with np.load(str(EXP/'energy_comparison/energy_predictions.npz')) as f: expected=f['u_anchored_energy']
        replay=float(np.max(np.abs(unew-expected)))
        report['max_center_replay_difference']=replay
        if replay>2e-6:raise RuntimeError('Frozen center prediction replay failed: '+str(replay))
        for k in ['label','source','u_true','truth_valid']:out[k]=a['eval_'+k]
    else:
        out['source']=a['source'];out['outer_holdout']=a['outer_holdout']
    print(kind,'center predictions complete; loading material map',flush=True)
    material=rt.MaterialMap()
    n=len(unew)
    # Bounded batches keep the native columns below 70 MB while retaining
    # enough parallel NumPy work for the shared login/analysis host.
    geometries=[];raw_parts=[];new_parts=[]
    for start in range(0,n,16000):
        sl=slice(start,min(start+16000,n))
        geo=material.build_columns(a['tk_pos'][sl],a['tk_dir'][sl])
        raw_parts.append(rt.propagate(d['u'][sl],d['R'][sl],geo))
        new_parts.append(rt.propagate(unew[sl],d['R'][sl],geo))
        quality=out['rich_quality'][sl]
        small={k:v[quality] if isinstance(v,np.ndarray) and k not in ['segment','faces'] else v for k,v in geo.items()}
        small['event_rows']=np.flatnonzero(quality)+start
        geometries.append(small)
        print(kind,'transport',sl.stop,'/',n,flush=True)
    raw=concatenate(raw_parts);new=concatenate(new_parts)
    for key in ['beta_planes','loss_GeV','valid','mass_GeV','initial_energy_GeV','final_energy_GeV']:
        out['raw_'+key]=raw[key];out['new_'+key]=new[key]
    del material
    good=raw['valid']&new['valid']
    out['common_transport_valid']=good
    report['transport_valid']=int(good.sum())
    report['unphysical_center_raw']=int((d['u']<=1).sum())
    report['unphysical_center_new']=int((unew<=1).sum())
    geometry={k:geometries[0][k] if k in ['segment','faces'] else np.concatenate([g[k] for g in geometries]) for k in geometries[0]}
    g_rows=geometry.pop('event_rows')
    prednaf=np.zeros(n,bool);prednaf[g_rows]=geometry['NaF_predicted']
    report['radiator_assignment_disagreement']=int(np.sum(out['rich_quality']&(prednaf!=out['NaF'])))
    report['target_z_difference_on_RICH_cm']=stats(geometry['target_z']-a['rich_pos'][g_rows,2])
    rich=out['rich_quality']&good
    out['rich_comparison_sample']=rich
    ref=out['rich_beta']
    br=raw['beta_planes'][:,-1];bn=new['beta_planes'][:,-1]
    report['rich_all']=compare(br[rich],bn[rich],ref[rich]) if rich.sum()>1 else {'n':int(rich.sum())}
    report['rich_naive_center_comparison']=compare(1/d['u'][rich],1/unew[rich],ref[rich]) if rich.sum()>1 else {}
    report['rich_bootstrap']=bootstrap(br[rich],bn[rich],ref[rich],out['source'][rich],out['label'][rich] if kind=='mc' else None)
    groups={
        'NaF':out['NaF'],'aerogel':~out['NaF'],'strict':out['rich_strict'],
        'beta_0.78_0.85':(ref>=.78)&(ref<.85),'beta_0.85_0.90':(ref>=.85)&(ref<.90),
        'beta_0.90_0.95':(ref>=.90)&(ref<.95),'beta_0.95_1.02':(ref>=.95)&(ref<=1.02)}
    if kind=='mc':groups.update({'C12':out['label']==0,'C13':out['label']==1})
    else:groups['outer_holdout']=out['outer_holdout']
    report['rich_groups']={}
    for name,mask in groups.items():
        sel=rich&mask
        if sel.sum()>1:
            report['rich_groups'][name]=compare(br[sel],bn[sel],ref[sel])
            if name in ['NaF','aerogel','C12','C13','outer_holdout','strict']:
                report['rich_groups'][name]['bootstrap']=bootstrap(br[sel],bn[sel],ref[sel],out['source'][sel],out['label'][sel] if kind=='mc' else None)
        else:report['rich_groups'][name]={'n':int(sel.sum())}
    report['expected_loss_RICH_sample']={name:stats(new['loss_GeV'][rich,j]) for j,name in enumerate(rt.LOSS_NAMES)}
    report['expected_total_loss_RICH_sample']=stats(new['loss_GeV'][rich].sum(axis=1))
    report['expected_beta_drop_RICH_sample']=stats(1/unew[rich]-bn[rich])
    sensitivity={}
    for name,kw in [('material_0.8',{'material_scale':.8}),('material_1.2',{'material_scale':1.2}),
                    ('common_mass_C12',{'mass_override':11.17793}),('common_mass_C13',{'mass_override':12.11221})]:
        pr=rt.propagate(d['u'][g_rows],d['R'][g_rows],geometry,**kw)
        pn=rt.propagate(unew[g_rows],d['R'][g_rows],geometry,**kw)
        mask=good[g_rows]&pr['valid']&pn['valid']
        if mask.sum()>1:sensitivity[name]=compare(pr['beta_planes'][mask,-1],pn['beta_planes'][mask,-1],ref[g_rows][mask])
    report['transport_sensitivity']=sensitivity
    if kind=='mc':
        tv=out['truth_valid']&np.isfinite(out['u_true'])&(out['u_true']>0)
        true=1/out['u_true']
        report['mc_center_all']=compare(1/d['u'][tv],1/unew[tv],true[tv])
        report['mc_center_bootstrap']=bootstrap(1/d['u'][tv],1/unew[tv],true[tv],out['source'][tv],out['label'][tv])
        sel=tv&rich
        report['mc_center_same_RICH_events']=compare(1/d['u'][sel],1/unew[sel],true[sel])
        report['mc_center_same_RICH_bootstrap']=bootstrap(1/d['u'][sel],1/unew[sel],true[sel],out['source'][sel],out['label'][sel])
        out['truth_beta_center']=true
        mom=a['eval_mevmom1'].astype(float);xyz=a['eval_mevcoo1'].astype(float)
        mass=np.where(out['label']==0,11.17793,12.11221)
        energy=np.sqrt(mom*mom+mass[:,None]**2)
        checkpoint=np.array([10,11,12,13,17])
        plane=np.array([0,1,2,3,6])
        realplanes=np.array([1.698,-2.318,-25.212,-29.228,-69.975])
        validtruth=(mom[:,checkpoint]>0)&np.isfinite(mom[:,checkpoint])&(np.abs(xyz[:,checkpoint,2]-realplanes[None,:])<.06)
        primary=(a['eval_mtrpri'][:,4]==1)&(a['eval_mtrmom'][:,4]>0)&(np.abs(a['eval_mtrcoo'][:,4,2]-1.698)<.03)
        report['mc_primary_L5_count']=int(primary.sum())
        beta_primary=a['eval_mtrmom'][:,4]/np.sqrt(a['eval_mtrmom'][:,4]**2+mass*mass)
        sel=primary&np.isfinite(beta_primary)&(beta_primary>0)
        if sel.sum()>1:report['mc_center_primary_L5']=compare(1/d['u'][sel],1/unew[sel],beta_primary[sel])
        out['mc_primary_L5']=primary
        report['energy_truth_closure']={}
        for j in range(1,len(checkpoint)):
            lo,hi=checkpoint[0],checkpoint[j]
            sel=good&validtruth[:,0]&validtruth[:,j]&primary
            truthloss=energy[:,lo]-energy[:,hi]
            predicted=new['loss_GeV'][:,:plane[j]].sum(axis=1)
            name='L5_to_'+rt.PLANE_NAMES[plane[j]]
            result={'truth_loss_GeV':stats(truthloss[sel]),'predicted_loss_GeV':stats(predicted[sel]),
                    'prediction_minus_truth_GeV':stats((predicted-truthloss)[sel]),
                    'correlation':float(np.corrcoef(truthloss[sel],predicted[sel])[0,1]),
                    'negative_truth_loss':int((truthloss[sel]<0).sum()),'truth_loss_over_1GeV':int((truthloss[sel]>1).sum())}
            truthbeta=mom[:,hi]/energy[:,hi]
            result['beta_at_checkpoint']=compare(raw['beta_planes'][sel,plane[j]],new['beta_planes'][sel,plane[j]],truthbeta[sel])
            report['energy_truth_closure'][name]=result
            out['truth_loss_'+name]=np.where(sel,truthloss,np.nan)
            if j==len(checkpoint)-1:
                out['truth_beta_post_S4']=np.where(sel,truthbeta,np.nan)
                result['RICH_sample']=stats((predicted-truthloss)[sel&rich])
                result['RICH_minus_postS4_truth_beta']=stats((ref-truthbeta)[sel&rich])
                # Substitution diagnosis: exact center beta/R only assess the
                # transport error. This is not a deployable model or training.
        paired=good&validtruth[:,0]&validtruth[:,-1]&primary&rich
        # Exact-center diagnostic uses frozen geometry only and true momentum
        # after all deployable predictions have been saved in memory.
        pg=paired[g_rows]
        truth_geo={k:(v[pg] if k not in ['segment','faces'] else v) for k,v in geometry.items()}
        truth_pr=rt.propagate(out['u_true'][g_rows][pg],mom[g_rows,10][pg]/6,truth_geo)
        if pg.sum()>1:
            report['transport_exact_center_diagnostic']={
                'predicted_RICH_minus_measured_RICH_beta':stats(truth_pr['beta_planes'][:,-1]-ref[g_rows][pg]),
                'predicted_postS4_minus_truth_beta':stats(truth_pr['beta_planes'][:,6]-mom[g_rows,17][pg]/energy[g_rows,17][pg])}
    np.savez_compressed(str(dest/(kind+'_predictions.npz')),**out)
    dump(dest/(kind+'_results.json'),report)
    rows=[]
    for group in ['rich_all','mc_center_all','mc_center_same_RICH_events']:
        if group not in report:continue
        for name in ['baseline_beta','new_beta','baseline_u','new_u']:
            if name in report[group]:rows.append({'group':group,'method':name,**report[group][name]})
    with (dest/(kind+'_metrics.tsv')).open('w') as f:
        if rows:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    print('DONE',kind,'selected',n,'RICH',int(rich.sum()),flush=True)
    print(json.dumps({k:report[k] for k in ['rich_all','rich_bootstrap','mc_center_all','mc_center_same_RICH_events'] if k in report},indent=2,default=lambda x:x.tolist()),flush=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('kind',choices=['mc','data']);ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--output',default='results/rich_transport');ap.add_argument('--cache',default='/tmp/betaEstimator-rich-transport')
    args=ap.parse_args();dest=ROOT/args.output;freeze(dest)
    cache=Path(args.cache);cache.mkdir(parents=True,exist_ok=True)
    path=cache/(args.kind+'_evaluation_input.npz')
    if path.exists():
        with np.load(str(path)) as f:a={k:f[k] for k in f.files}
        print('Loaded completed evaluation cache',path,len(a['run']),flush=True)
    else:
        a=extract(args.kind,args.workers,dest)
        np.savez_compressed(str(path),**a)
    evaluate(args.kind,a,dest)

if __name__=='__main__':main()
