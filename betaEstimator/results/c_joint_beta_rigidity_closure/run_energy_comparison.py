import sys, json, math, csv, time, multiprocessing as mp
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import uproot
HERE=Path(__file__).resolve().parent
SAVE=HERE/'expanded_validation'/'energy_comparison'
BASE=HERE.parent/'c_isotope_proportion_closure'; OUT=HERE/'expanded_validation'
sys.path[:0]=[str(BASE),str(HERE)]
import run_closure as rc
import run_energy_residualized as er
import run_anchored_joint_scan as aj
import run_joint_identifiability_scan as js
import plot_beta_residual_trends as trends
C_LIGHT=rc.C_LIGHT

def one_file(arg):
    label, source, path, rows, target_run, target_event = arg
    names=['run','event','tof_tl','tof_leng','tof_pos','tof_barid','tof_edep','tk_ql']
    with uproot.open(path) as f:
        a=f['amstreea'].arrays(names,library='np')
    ids=np.rec.fromarrays([a['run'].astype(np.int64),a['event'].astype(np.int64)],names='run,event')
    tids=np.rec.fromarrays([target_run.astype(np.int64),target_event.astype(np.int64)],names='run,event')
    common, ie, it=np.intersect1d(ids,tids,assume_unique=False,return_indices=True)
    if len(common)!=len(rows):
        raise RuntimeError('match failure label=%d source=%d got=%d need=%d'%(label,source,len(common),len(rows)))
    order=np.argsort(it); ie=ie[order]; rows_out=rows[it[order]]
    times=a['tof_tl'][ie].astype(np.float64)
    path_ns=a['tof_leng'][ie].astype(np.float64)/C_LIGHT
    cp=path_ns-path_ns.mean(axis=1,keepdims=True)
    ct=times-times.mean(axis=1,keepdims=True)
    den=np.sum(cp*cp,axis=1)
    u=np.sum(cp*ct,axis=1)/den
    res=ct-u[:,None]*cp
    span=np.ptp(path_ns,axis=1)
    pos=a['tof_pos'][ie].astype(np.float64)
    dz=pos[:,3,2]-pos[:,0,2]
    dx=(pos[:,3,0]-pos[:,0,0])/dz; dy=(pos[:,3,1]-pos[:,0,1])/dz
    continuous=np.column_stack([
        pos[:,:,0]/120.0,pos[:,:,1]/120.0,dx,dy,cp,res/span[:,None],
        np.log(a['tof_edep'][ie].astype(np.float64)/(rc.CHARGE*rc.CHARGE)),
        2.0*np.log(a['tk_ql'][ie,1:8].astype(np.float64)/rc.CHARGE),
        np.log(span/4.0)]).astype(np.float64)
    return rows_out, continuous, a['tof_barid'][ie].astype(np.int8), u

def main():
    started=time.time()
    # Train/freeze candidate on original 7200-event weak split only.
    train,valid=er.load_cache(BASE/'controlled_cache.npz')
    base_state=rc.make_feature_state(train,'no_edep'); bt=rc.apply_feature_state(train,base_state)
    response=er.response_state(train); raw_t=er.energy_raw(train,'fitted',response)
    std=er.standardizer(raw_t); zt=er.apply_standardizer(raw_t,std)
    dt=np.column_stack([zt,zt*zt]); dcenter=np.mean(dt,axis=0); dscale=np.std(dt,axis=0); dscale[dscale<1e-6]=1
    dt=np.clip((dt-dcenter)/dscale,-5,5); xt=np.column_stack([bt,dt])
    bv=rc.apply_feature_state(valid,base_state); raw_v=er.energy_raw(valid,'fitted',response); zv=er.apply_standardizer(raw_v,std)
    dv=np.column_stack([zv,zv*zv]); dv=np.clip((dv-dcenter)/dscale,-5,5); xv=np.column_stack([bv,dv])
    xa_t,_,astate=aj.anchored_projection(xt[:,1:],xv[:,1:],train,valid,'beta_fitted')
    weak={k:train[k] for k in ['u','R','relR','bin']}; density=rc.MixtureDensity(weak,xa_t,train['prior_c13'])
    logs=np.asarray(json.loads((BASE/'weak_frozen_models.json').read_text())['models']['raw']['logsigma'])
    p=xa_t.shape[1]; init=np.r_[np.zeros(p),logs]
    bounds=[(0,0)]+[(-.05,.05)]*(p-1)+[(math.log(.0015),math.log(.08))]*rc.N_BINS
    fit=minimize(density.value_gradient,init,jac=True,method='L-BFGS-B',bounds=bounds,
                 options={'maxiter':5000,'maxls':100,'ftol':1e-14,'gtol':5e-7,'maxcor':30})
    assert fit.success, fit.message
    (SAVE/'energy_frozen_model.json').write_text(json.dumps(aj.jsonable({'base_state':base_state,'response':response,'standardizer':std,'derived_center':dcenter,'derived_scale':dscale,'anchor':astate,'parameters':fit.x}),indent=2)+'\n')
    theta=fit.x[:p]; grad=float(np.max(np.abs(density.value_gradient(fit.x)[1])))
    # Existing independent expanded event list and frozen comparisons.
    with np.load(OUT/'predictions.npz',allow_pickle=False) as f: old={k:f[k] for k in f.files}
    n=len(old['u']); continuous=np.empty((n,30),dtype=np.float64); bars=np.empty((n,4),dtype=np.int8); ucheck=np.empty(n)
    manifest={0:{},1:{}}; lm={'C12':0,'C13':1}
    with (OUT/'manifest.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'): manifest[lm[r['isotope']]][int(r['source'])]=r['path']
    tasks=[]
    for label in (0,1):
      for source,path in sorted(manifest[label].items()):
        rows=np.flatnonzero((old['label']==label)&(old['source']==source))
        tasks.append((label,source,path,rows,old['run'][rows],old['event'][rows]))
    with mp.Pool(processes=6) as pool:
      for j,(rows,c,b,u) in enumerate(pool.imap_unordered(one_file,tasks,chunksize=1),1):
        continuous[rows]=c; bars[rows]=b; ucheck[rows]=u
        if j%48==0: print('files',j,'/',len(tasks),'elapsed',time.time()-started,flush=True)
    maxdu=float(np.max(np.abs(ucheck-old['u'])))
    if maxdu>2e-6: raise RuntimeError('u replay '+str(maxdu))
    d={'u':old['u'],'R':old['R'],'bin':old['bin'],'continuous':continuous,'bars':bars}
    b=rc.apply_feature_state(d,base_state); raw_e=er.energy_raw(d,'fitted',response); z=er.apply_standardizer(raw_e,std)
    der=np.column_stack([z,z*z]); der=np.clip((der-dcenter)/dscale,-5,5)
    x=np.column_stack([b,der])[:,1:]
    proj=(x-astate['bin_centers'][d['bin']]).dot(astate['null_projection'])/astate['scale']
    xa=np.column_stack([np.zeros(len(proj)),proj]); corr=xa.dot(theta); u_new=old['u']-corr
    truth=old['truth_valid'] & np.isfinite(old['u_true']) & (old['u_true']>0)
    methods={'raw':old['u_raw'],'beta_only_frozen':old['u_beta_only_frozen'],'joint_anchored':old['u_joint_anchored'],'anchored_energy':u_new}
    np.savez_compressed(str(SAVE/'energy_predictions.npz'),run=old['run'],event=old['event'],label=old['label'],source=old['source'],u_anchored_energy=u_new)
    summary={}
    for name,u in methods.items():
      du=u[truth]-old['u_true'][truth]; db=1/u[truth]-1/old['u_true'][truth]
      m2=(6*old['R'][truth])**2*(u[truth]**2-1); labels=old['label'][truth]; m2err=m2-rc.MASSES[labels]**2
      summary[name]={'median_du':float(np.median(du)),'w68_du':rc.w68(du),'median_dbeta':float(np.median(db)),
                     'w68_dbeta':rc.w68(db),'mass2_auc':rc.auc(m2,labels),'rejection50':js.roc_rejection(m2,labels),
                     'median_mass2_error':float(np.median(m2err)),'w68_mass2_error':rc.w68(m2err)}
    rw=summary['raw']['w68_du']
    for v in summary.values(): v['du_gain_vs_raw_pct']=100*(1-v['w68_du']/rw)
    detail=[]
    for bid in [-1,0,1,2,3]:
      for lab in [-1,0,1]:
        sel=truth & ((old['bin']==bid) if bid>=0 else True) & ((old['label']==lab) if lab>=0 else True)
        du=u_new[sel]-old['u_true'][sel]
        detail.append({'bin':bid,'label':lab,'n':int(sel.sum()),'median_du':float(np.median(du)),'w68_du':rc.w68(du)})
    beta_truth=1/old['u_true'][truth]; trend=[]
    for name,u in methods.items():
      rr=1000*(1/u[truth]-beta_truth)
      for lo,hi in zip(trends.EDGES[:-1],trends.EDGES[1:]):
        vals=rr[(beta_truth>=lo)&(beta_truth<hi)]; f=trends.fit_core(vals)
        trend.append({'method':name,'beta_center':float((lo+hi)/2),'n':len(vals),**f})
    eraw=old['u'][truth]-old['u_true'][truth]; enew=u_new[truth]-old['u_true'][truth]; lab=old['label'][truth]; src=old['source'][truth]
    group={(l,s):np.flatnonzero((lab==l)&(src==s)) for l in (0,1) for s in np.unique(src[lab==l])}
    keys={l:sorted(k for k in group if k[0]==l) for l in (0,1)}; rng=np.random.RandomState(290929); gains=[]
    for _ in range(300):
      idx=np.concatenate([group[keys[l][i]] for l in (0,1) for i in rng.randint(0,len(keys[l]),size=len(keys[l]))])
      gains.append(100*(1-rc.w68(enew[idx])/rc.w68(eraw[idx])))
    result={'training':{'n':len(train['u']),'features':int(p),'success':bool(fit.success),'iterations':int(fit.nit),'objective':float(fit.fun),'gradient_inf':grad},
            'expanded':{'selected':n,'truth_valid':int(truth.sum()),'files':len(tasks),'u_replay_max':maxdu,
                        'correction_rms':rc.rms(corr),'correction_median':float(np.median(corr)),'correction_label_auc':rc.auc(corr,old['label']),
                        'summary':summary,'new_detail':detail,'beta_trend':trend,
                        'bootstrap_gain_pct':{'median':float(np.median(gains)),'lo95':float(np.percentile(gains,2.5)),'hi95':float(np.percentile(gains,97.5))}},
            'elapsed_seconds':time.time()-started}
    out=SAVE/'probe_results.json'; out.write_text(json.dumps(result,indent=2)+'\n')
    print('DONE',out,'elapsed',result['elapsed_seconds'],'gain',summary['anchored_energy']['du_gain_vs_raw_pct'],result['expanded']['bootstrap_gain_pct'],flush=True)
if __name__=='__main__': main()
