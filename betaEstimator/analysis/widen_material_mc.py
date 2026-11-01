#!/usr/bin/env python3
import argparse
import csv
import json
import io
import os
from pathlib import Path
import subprocess
import time
import numpy as np
import uproot
import rich_transport as rt
import evaluate_rich_transport as ev
import flight_beta_wide as wide

OUT=rt.ROOT/'results/time_local_ratio/material_closure_wide_R'
MANIFEST=ev.EXP/'manifest.tsv'
MASS=[11.17793,12.11221]


def prepare():
    OUT.mkdir(exist_ok=True)
    for name in ['parts','condor']: (OUT/name).mkdir(exist_ok=True)
    records=list(csv.DictReader(MANIFEST.open(),delimiter='\t'))
    ev.dump(OUT/'protocol.json',dict(records=records,R_GV=[2.5,200.],
        reference='Same loose RICH quality plus finite 4 < rich_q[1] < 8',
        transport='Same fixed nominal C12 mean inverse-beta increment from L5 to RICH',
        truth='mevmom1[:,10] with the production isotope mass, same as the original MC evaluation',
        unchanged_flight_model=True,MC_only=True,
        selection_scope='Detector and RICH quality, fixed nominal material validity; no learned/raw TOF forward-transport-validity gate, since this check uses no TOF beta prediction',
        manifest_sha256=ev.digest(MANIFEST),
        script_sha256=ev.digest(__file__)))
    ads=json.loads(subprocess.check_output(['condor_status','-json']))
    eligible=sorted(set(a['Machine'] for a in ads if a.get('OpSysAndVer')=='CentOS7' and a['Machine']!='amswn135.ams.sdiat.cn'))
    free=sorted(set(a['Machine'] for a in ads if a['Machine'] in eligible and a.get('State')=='Unclaimed' and a.get('Activity')=='Idle'))
    hosts=free+[v for v in eligible if v not in free]
    lines=[]
    for worker in range(50):
        lines.append('%d %s'%(worker,';'.join(hosts[worker::50])))
    script=OUT/'condor/run.sh'
    script.write_text('#!/bin/bash\nset -eo pipefail\ncd /home/ams/hxwu/AMSSoft/betaEstimator\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\npython3 analysis/widen_material_mc.py worker --worker "$1"\n')
    sub=OUT/'condor/jobs.sub'
    sub.write_text('universe = vanilla\nexecutable = /bin/bash\narguments = '+str(script)+' $(worker)\ngetenv = true\nrequest_cpus = 32\nrequest_memory = 8GB\nshould_transfer_files = NO\nrequirements = (TARGET.OpSysAndVer == "CentOS7") && stringListMember(TARGET.Machine, "$(allowed_hosts)", ";")\noutput = '+str(OUT)+'/condor/$(Cluster)_$(Process).out\nerror = '+str(OUT)+'/condor/$(Cluster)_$(Process).err\nlog = '+str(OUT)+'/condor/$(Cluster)_$(Process).log\nqueue worker,allowed_hosts from (\n'+'\n'.join(lines)+'\n)\n')
    text=subprocess.check_output(['condor_submit','-name','scheduler@schedd01.ams.sdiat.cn',str(sub)],universal_newlines=True)
    import re
    cluster=int(re.search(r'cluster (\d+)',text).group(1));ev.dump(OUT/'submission.json',dict(cluster=cluster,workers=50,eligible_hosts=eligible));print(text,flush=True)


def worker(index):
    if not hasattr(np, 'isin'):
        np.isin = lambda values, choices: np.in1d(values, choices).reshape(np.shape(values))
    records=json.loads((OUT/'protocol.json').read_text())['records']
    wide.set_edges(wide.R_EDGES)
    material=rt.MaterialMap()
    for j in range(index,len(records),50):
        target=OUT/'parts'/('%04d.npz'%j)
        if target.exists():continue
        record=records[j];label=0 if record['isotope']=='C12' else 1
        branches=list(dict.fromkeys(ev.MEASURED+ev.SELECTION+['rich_q','mevmom1','mpar','mch']))
        with uproot.open(record['path']) as source:
            tree=source['amstreea'];a=tree.arrays(branches,library='np');raw=tree.num_entries
        d=rt.center_features(a);R=d['R'];q=a['rich_q'][:,1]
        mask=wide.selection(a)&(a['mch']==6)&(a['mpar']==(67 if label==0 else 117))
        # Match the original detector-only geometry and timing validity checks.
        dz=a['tof_pos'][:,3,2]-a['tof_pos'][:,0,2]
        mask &= np.isfinite(a['tof_pos']).all((1,2)) & (np.abs(dz)>1.)
        mask &= (np.ptp(a['tof_leng'],axis=1)>0)
        mask &= rt.rich_quality(a)&np.isfinite(q)&(q>4)&(q<8)
        rows=np.flatnonzero(mask);b=ev.sub(a,rows);r=R[rows]
        geo=material.build_columns(b['tk_pos'],b['tk_dir'])
        nominal=np.sqrt(1+(MASS[0]/(6*r))**2)
        prop=rt.propagate(nominal,r,geo,mass_override=MASS[0])
        delay=1/prop['beta_planes'][:,-1]-nominal
        mom=b['mevmom1'][:,10].astype(float)
        refu=1/b['rich_beta'][:,0].astype(float)-delay
        good=prop['valid']&np.isfinite(mom)&(mom>0)&np.isfinite(refu)&(refu>0)
        out=dict(run=b['run'][good],event=b['event'][good],R=r[good],
            beta_mc=mom[good]/np.hypot(mom[good],MASS[label]),u_reference=refu[good],
            isotope=np.full(int(good.sum()),label,np.int8),NaF=b['rich_NaF'][good].astype(bool),
            source=np.full(int(good.sum()),int(record['source']),np.int32),delay_u=delay[good])
        archive=io.BytesIO()
        np.savez_compressed(archive,**out)
        target.with_suffix('.partial.npz').write_bytes(archive.getvalue())
        target.with_suffix('.partial.npz').replace(target)
        ev.dump(target.with_suffix('.json'),dict(record=record,raw=raw,reference_candidates=len(rows),n=int(good.sum())))
        print('PART',j,'N',int(good.sum()),'NaF',int(out['NaF'].sum()),flush=True)


def assemble():
    records=json.loads((OUT/'protocol.json').read_text())['records']
    paths=[OUT/'parts'/('%04d.npz'%j) for j in range(len(records))]
    if not all(p.exists() for p in paths):raise RuntimeError('Incomplete wide R extraction')
    chunks=[]
    for path in paths:
        with np.load(str(path)) as f: chunks.append({k:f[k] for k in f.files})
    a=ev.concatenate(chunks)
    ids=np.rec.fromarrays([a['isotope'],a['source'],a['run'],a['event']])
    if len(np.unique(ids))!=len(ids):raise RuntimeError('Duplicate MC identities')
    with uproot.recreate(str(OUT/'comparisons.root')) as root:
        for naf,rad in [(True,'NaF'),(False,'AGL')]:root['mc_'+rad]=ev.sub(a,a['NaF']==naf)
    stats={}
    for naf,rad in [(True,'NaF'),(False,'AGL')]:
        s=a['NaF']==naf;res=1/a['u_reference'][s]-a['beta_mc'][s]
        stats[rad]=dict(n=int(s.sum()),R_range=[float(a['R'][s].min()),float(a['R'][s].max())],residual=ev.stats(res))
    # Verify that widening R preserves every previously evaluated event.
    old=uproot.open(str(OUT.parent/'material_closure/comparisons.root'));overlap={}
    for naf,rad in [(True,'NaF'),(False,'AGL')]:
        b=old['mc_'+rad].arrays(library='np');sel=(a['NaF']==naf)&(a['R']<6)
        oldids=np.rec.fromarrays([b[k] for k in ['isotope','source','run','R']])
        newids=np.rec.fromarrays([a[k][sel] for k in ['isotope','source','run','R']])
        if len(np.unique(oldids))!=len(oldids) or len(np.unique(newids))!=len(newids):raise RuntimeError('Ambiguous old-cache association')
        common,i,j=np.intersect1d(oldids,newids,return_indices=True)
        if len(common)!=len(oldids):raise RuntimeError('Original low-R events missing '+rad+' '+str((len(oldids),len(newids),len(common))))
        delta=np.max(np.abs(b['u_reference'][i]-a['u_reference'][sel][j]))
        if delta>1e-7:raise RuntimeError('Low-R transport replay differs')
        overlap[rad]=dict(old_n=len(oldids),wide_low_R_n=len(newids),additional_without_TOF_forward_transport_gate=len(newids)-len(oldids),join_fields=['isotope','source','run','R'],unique_join_verified=True,max_u_reference_difference=float(delta))
    ev.dump(OUT/'verification.json',dict(passed=True,R_GV=[2.5,200.],stats=stats,low_R_replay=overlap,files=len(paths)))
    print('WIDE_R_VERIFIED',stats,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','worker','assemble']);p.add_argument('--worker',type=int);a=p.parse_args()
    if a.stage=='prepare':prepare()
    elif a.stage=='worker':worker(a.worker)
    else:assemble()
