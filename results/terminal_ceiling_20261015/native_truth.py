#!/usr/bin/env python3
import argparse,json,subprocess,time,hashlib
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np,uproot
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
REF=ROOT/'results/flight_calibration_research/conditional_mc/reference'
ENV=ROOT.parent/'AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh'
RAW=Path('/lustre03/data/eos/ams/MC/AMS02/2022/C.B1308/c12.pl1.l1.612000.6_05')

def prepare():
 a=np.load(str(HERE/'broad_Z6.npz'));manifest=json.loads((HERE/'broad_Z6_manifest.json').read_text());out=HERE/'native_parts';out.mkdir(exist_ok=True);records=[]
 for r in manifest['sources']:
  sid=r['source'];select=a['source']==sid;path=Path(r['path']);raw=RAW/(path.stem.replace('_','.')+'.root')
  if not raw.exists():raise FileNotFoundError(str(raw))
  b=uproot.open(str(path))['amstreea'].arrays(['run','event','itrtrack','ibetah'],library='np')
  key=b['run'].astype(np.uint64)*np.uint64(2**32)+b['event'].astype(np.uint64);order=np.argsort(key)
  ids=a['run'][select].astype(np.uint64)*np.uint64(2**32)+a['event'][select].astype(np.uint64);j=np.searchsorted(key[order],ids)
  if not np.array_equal(key[order[j]],ids):raise ValueError('Native target identity mismatch')
  ix=order[j];target=out/('%03d.txt'%sid)
  np.savetxt(str(target),np.column_stack([a['run'][select],a['event'][select],b['itrtrack'][ix],b['ibetah'][ix],a['tof_barid'][select]]),fmt='%d')
  records.append(dict(source=sid,path=str(raw),target=str(target),output=str(out/('%03d.root'%sid)),events=int(select.sum()),role=r['role']))
 plan=dict(records=records,events=sum(r['events'] for r in records),decoder=str(REF/'extract_parent'),decoder_source=str(ROOT/'analysis/conditional_mc_reference_parent.C'),decoder_sha256=hashlib.sha256((REF/'extract_parent').read_bytes()).hexdigest(),selection='Reuse complete measured+primary-L5 broad C12 cohort. Preserve all truth records; primary TOF selection occurs by track/particle/parent identity, not error magnitude.')
 (HERE/'native_plan.json').write_text(json.dumps(plan,indent=2));print('NATIVE_PLAN',len(records),plan['events'],flush=True)

def one(r):
 output=Path(r['output']);note=output.with_suffix('.json')
 if note.exists() and json.loads(note.read_text()).get('complete'):return json.loads(note.read_text())
 if output.exists():raise RuntimeError('Unverified partial native output '+str(output))
 command=['bash','-c','source "$1" >/dev/null; exec "$2" "$3" "$4" "$5"','_',str(ENV),str(REF/'extract_parent'),r['path'],r['target'],r['output']]
 start=time.monotonic()
 with output.with_suffix('.log').open('w') as log:proc=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
 result=dict(source=r['source'],events=r['events'],complete=proc.returncode==0,exit_code=proc.returncode,seconds=time.monotonic()-start)
 if proc.returncode==0:
  b=uproot.open(r['output'])['reference'].arrays(['run','event','found','track_exists','beta_exists','beta_track','requested_track'],library='np')
  expected=np.loadtxt(r['target'],dtype=np.int64,ndmin=2)
  result['identities_equal']=bool(np.array_equal(b['run'],expected[:,0]) and np.array_equal(b['event'],expected[:,1]))
  result['all_found']=bool(b['found'].all());result['associations_valid']=int(np.sum(b['track_exists']&b['beta_exists']&(b['beta_track']==b['requested_track'])))
  result['complete']&=result['identities_equal'] and result['all_found']
 note.write_text(json.dumps(result,indent=2));print('NATIVE_SOURCE',result,flush=True)
 if not result['complete']:raise RuntimeError(str(result))
 return result

def run(workers):
 plan=json.loads((HERE/'native_plan.json').read_text())
 with ThreadPoolExecutor(max_workers=workers) as pool:results=list(pool.map(one,plan['records']))
 (HERE/'native_complete.json').write_text(json.dumps(dict(complete=True,sources=results),indent=2))
 print('NATIVE_COMPLETE',sum(r['events'] for r in results),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run']);p.add_argument('--workers',type=int,default=4);args=p.parse_args()
 if args.stage=='prepare':prepare()
 else:run(args.workers)
