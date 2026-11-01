#!/usr/bin/env python3
"""Prepare, run and verify an isolated native measured-readout supplement."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import uproot

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/flight_calibration_research/readout_audit'
FIELDS = ['run','event','itrtrack','ibetah','tof_barid','tof_adca','tof_edep',
          'tof_rtime','tof_rtime_ht','tof_rtime_sht','tof_pos','tk_pos','tk_dir',
          'tof_tl','tof_leng']
ENV = '/nas01/home/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')

def stats(a):
    a=np.asarray(a).ravel();a=a[np.isfinite(a)]
    if not len(a):return dict(n=0)
    q=np.percentile(a,[0,16,50,84,100])
    return dict(n=len(a),mean=float(a.mean()),q0=float(q[0]),q50=float(q[2]),q100=float(q[4]),w68=float((q[3]-q[1])/2))

def prepare(count):
    OUT.mkdir(parents=True,exist_ok=True)
    inventory=json.loads((ROOT/'results/flight_calibration_research/input_inventory.json').read_text())
    records=[]
    for rec in inventory['records']:
        source=int(rec['source']);part=ROOT/'results/time_local_beta/parts'/('%05d'%source)
        with np.load(str(part/'measured.npz'),allow_pickle=False) as f:
            keys=list(zip(f['run'].astype(int),f['event'].astype(int)))
        wanted=set(keys)
        with uproot.open(rec['path']) as f:
            tree=f['amstreea'];a=tree.arrays(FIELDS,library='np')
        ix=np.array([i for i,k in enumerate(zip(a['run'].astype(int),a['event'].astype(int))) if k in wanted],int)
        if not len(ix):raise ValueError('No native training targets')
        ix=ix[np.unique(np.linspace(0,len(ix)-1,min(count,len(ix))).astype(int))]
        dest=OUT/('%05d'%source);dest.mkdir(exist_ok=True)
        data={k:v[ix] for k,v in a.items()};np.savez_compressed(str(dest/'original.npz'),**data)
        targets=['%d %d %d %d'%v for v in zip(data['run'],data['event'],data['itrtrack'],data['ibetah'])]
        (dest/'targets.txt').write_text('\n'.join(targets)+'\n')
        runs=set(int(x) for x in data['run'])
        # The time-local list contains only newly supplemented events; the
        # validation list supplies native inputs for reused measured events.
        lists=[part/'raw_files.txt',ROOT/'results/tofcalib_validation/parts'/('%05d'%source)/'raw_files.txt']
        available=sorted(set(p for name in lists for p in name.read_text().splitlines() if p.strip()))
        paths=[p for p in available if int(Path(p).name.split('.')[0]) in runs]
        missing=runs-set(int(Path(p).name.split('.')[0]) for p in paths)
        if missing:raise ValueError('Missing native runs '+str(sorted(missing)))
        (dest/'raw_files.txt').write_text('\n'.join(paths)+'\n')
        records.append(dict(source=source,path=rec['path'],events=len(ix),runs=sorted(runs),native_files=len(paths),targets_sha256=sha(dest/'targets.txt')))
        print('PREPARED',source,len(ix),'events',len(runs),'runs',flush=True)
    write(OUT/'protocol.json',dict(fields=FIELDS,records=records,selection='Uniform positions in existing detector-only training events per source; no new beta/RICH selection',
        RICH_read=False,MC_truth_read=False,raw_ADC_definition='TrClusterR::GetSignal(j,0); no optional signal corrections; stored cluster is already a reconstructed detector object',
        calibration_diagnostics='Tracker hit GetEdepH side 0/1/2 converted keV to MeV; TOF GetEdep dynode/anode/default; these are not declared raw or newly calibrated energy',
        script_sha256=sha(__file__),native_source_sha256=sha(ROOT/'analysis/flight_readout_supplement.C'),binary_sha256=sha(OUT/'supplement')))

def run_one(rec):
    dest=OUT/('%05d'%rec['source']);result=dest/'native.root'
    if result.exists():return dict(source=rec['source'],reused=True)
    partial=dest/'native.partial.root'
    if partial.exists():raise RuntimeError('Inspect partial native output before retrying: '+str(partial))
    cmd='source %s >/dev/null; exec %s %s %s %s'%(ENV,OUT/'supplement',dest/'raw_files.txt',dest/'targets.txt',partial)
    with (dest/'native.log').open('w') as log:
        code=subprocess.call(['bash','-lc',cmd],cwd=str(ROOT),stdout=log,stderr=subprocess.STDOUT)
    note=dict(source=rec['source'],exit_code=code)
    write(dest/'run.json',note)
    if code!=0:raise RuntimeError('Native extraction failed: '+str(note))
    partial.rename(result)
    print('NATIVE_COMPLETE',rec['source'],flush=True)
    return note

def verify():
    protocol=json.loads((OUT/'protocol.json').read_text());rows=[]
    all_run=[];all_event=[]
    for rec in protocol['records']:
        dest=OUT/('%05d'%rec['source'])
        with np.load(str(dest/'original.npz'),allow_pickle=False) as f:a={k:f[k] for k in f.files}
        with uproot.open(str(dest/'native.root')) as f:
            b=f['readout'].arrays(library='np')
        for key in ['run','event','itrtrack','ibetah']:
            if not np.array_equal(a[key],b[key]):raise ValueError('Association mismatch '+key)
        if not (np.all(b['selected']==1) and np.all(b['found']==1)):raise ValueError('Missing native events')
        comparison={}
        for key in ['tof_barid','tof_adca','tof_rtime','tof_rtime_ht','tof_rtime_sht','tof_edep']:
            orig=a[key];new=b[key] if key!='tof_edep' else b[key][:,:,2]
            finite=np.isfinite(orig)&np.isfinite(new)
            positive=finite&(orig>0)&(new>0)
            comparison[key]=dict(total=int(orig.size),equal=int((orig==new).sum()),max_abs=float(np.max(np.abs(orig[finite]-new[finite]))) if finite.any() else None,
                positive_n=int(positive.sum()),positive_max_abs=float(np.max(np.abs(orig[positive]-new[positive]))) if positive.any() else None)
        mismatch=[];max_strip_error=0.;nstrip=0;overflow=0;badstrip=0
        for e in range(len(a['run'])):
            for l in range(9):
                for s in range(2):
                    off=int(b['tk_offset'][e,l,s]);n=int(b['tk_nstrip'][e,l,s])
                    if off<0:continue
                    sig=np.asarray(b['strip_signal'][e][off:off+n],float)
                    diff=abs(sig[sig>0].sum()-b['tk_sum'][e,l,s])
                    max_strip_error=max(max_strip_error,float(diff))
                    if diff>1.e-3*max(1.,abs(float(b['tk_sum'][e,l,s]))):mismatch.append([e,l,s,diff])
                    nstrip+=n;overflow+=int(np.sum(b['strip_overflow'][e][off:off+n]));badstrip+=int(np.sum(np.asarray(b['strip_status'][e][off:off+n])!=0))
        if mismatch:raise ValueError('Raw ADC sum mismatch '+str(mismatch[:2]))
        dep=b['tof_edep'];default=dep[:,:,2];dynode=dep[:,:,0];anode=dep[:,:,1]
        good=np.isfinite(default)&(default>0)
        row=dict(source=rec['source'],n=len(a['run']),comparison=comparison,nstrip=nstrip,overflow_strips=overflow,bad_status_strips=badstrip,
            raw_positive_sum_max_abs_difference=max_strip_error,tof_default_equals_dynode=int((good&(default==dynode)).sum()),
            tof_default_equals_anode=int((good&(default==anode)).sum()),tof_positive_hits=int(good.sum()),
            tk_raw_sum=stats(b['tk_sum']),tk_calibrated_energy_MeV=stats(b['tk_edep']),native_sha256=sha(dest/'native.root'))
        rows.append(row);all_run.extend(a['run']);all_event.extend(a['event'])
        print('VERIFIED',rec['source'],comparison,flush=True)
    if len(set(zip(all_run,all_event)))!=len(all_run):raise ValueError('Duplicate identity')
    result=dict(complete=True,events=len(all_run),runs=len(np.unique(all_run)),association_mismatches=0,raw_strip_sum_check_passed=True,rows=rows,
        beta_resolution_test=False,absolute_energy_calibration=False,script_sha256=sha(__file__))
    write(OUT/'verification.json',result)
    print('AUDIT_COMPLETE',len(all_run),'events',len(np.unique(all_run)),'runs',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','verify','all']);p.add_argument('--events-per-source',type=int,default=64);p.add_argument('--workers',type=int,default=3);a=p.parse_args()
    if a.stage in ['prepare','all']:prepare(a.events_per_source)
    if a.stage in ['run','all']:
        records=json.loads((OUT/'protocol.json').read_text())['records']
        with ThreadPoolExecutor(max_workers=a.workers) as pool:list(pool.map(run_one,records))
    if a.stage in ['verify','all']:verify()
