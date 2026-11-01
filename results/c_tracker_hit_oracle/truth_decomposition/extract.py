#!/usr/bin/env python3
"""Identity-preserving extraction for truth / measurement decomposition. No fitting."""
import os,csv,json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import uproot
HERE=Path(__file__).resolve().parent
PARENT=HERE.parent
BRANCHES=['run','event','mmom','mch','mmass','mpar','mtrmom','mtrpar','mtrcoo','mtrcoo1','mtrpri','mtrz','mevcoo1','mevdir1','mevmom1','tk_hitc','tk_pos','tk_dir','tk_hitb','tk_lid','tk_cad','tk_res','tk_res1','tk_rigidity1','tk_rigidityk','tk_rigidityi','tk_erigidity1','tk_chis1','tk_adc','tk_cof','tk_edep','tk_qlqs','tof_edep','tof_leng','tof_tl','tof_betah','tof_betahmc']
def one(task):
    path,idx,run,evt=task
    with uproot.open(path) as f:
        t=f['amstreea']; names=[b for b in BRANCHES if b in t.keys()]
        ids=t.arrays(['run','event'],library='np')
        lookup={(int(r),int(e)):j for j,(r,e) in enumerate(zip(ids['run'],ids['event']))}
        loc=np.array([lookup[(int(r),int(e))] for r,e in zip(run,evt)])
        if len(lookup)!=len(ids['run']): raise RuntimeError('duplicate identities '+path)
        a=t.arrays(names,library='np')
        return idx,{k:v[loc] for k,v in a.items()}
def main():
    start=time.time(); z=np.load(str(PARENT.parent/'c_mass_toy_oracle/predictions.npz'))
    sel=z['label']==0
    data={k:z[k][sel] for k in ['run','event','source','label','bin','R','R_true','u','u_true','role']}
    data['event_index']=np.flatnonzero(sel).astype('int32')
    jobs=list(csv.DictReader(open(str(PARENT/'tkadc_production/jobs.tsv')),delimiter='\t'))
    tasks=[]
    for row in jobs:
        idx=np.flatnonzero(data['source']==int(row['source']))
        tasks.append((row['output_root'],idx,data['run'][idx],data['event'][idx]))
    filled=np.zeros(len(data['R']),bool); arrays={}
    with ProcessPoolExecutor(max_workers=4) as pool:
        for j,(idx,a) in enumerate(pool.map(one,tasks)):
            for k,v in a.items():
                if k not in arrays: arrays[k]=np.empty((len(filled),)+v.shape[1:],dtype=v.dtype)
                arrays[k][idx]=v
            filled[idx]=True
            if j%32==31: print('EXTRACT',j+1,len(tasks),'events',sum(filled),'seconds',round(time.time()-start,1),flush=True)
    assert filled.all()
    assert np.array_equal(arrays['run'],data['run']) and np.array_equal(arrays['event'],data['event'])
    assert np.array_equal(arrays['tk_rigidityk'][:,1,1],data['R'])
    for k in data: arrays['ref_'+k]=data[k]
    np.savez_compressed(str(HERE/'aligned.npz'),**arrays)
    # Old ROOT can consume this explicit TTree; no RNTuple dependency.
    payload={k:v for k,v in arrays.items() if k not in ('ref_role',)}
    payload['role']=np.where(data['role']==b'train',0,np.where(data['role']==b'valid',1,2)).astype('int32')
    with uproot.recreate(str(HERE/'aligned.root')) as f:
        schema={k:(v.dtype if v.ndim==1 else np.dtype((v.dtype,v.shape[1:]))) for k,v in payload.items()}
        f.mktree('events',schema);f['events'].extend(payload)
    summary={'events':len(filled),'sources':len(jobs),'roles':{s:int(np.sum(data['role']==s.encode())) for s in ['train','valid','test']},'branches':list(arrays),'replay':'exact','seconds':time.time()-start}
    (HERE/'extraction.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
