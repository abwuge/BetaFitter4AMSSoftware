#!/usr/bin/env python3
"""Continue the authorized local-time experiment through final verification."""
import json
import multiprocessing as mp
from pathlib import Path
import re
import subprocess
import sys
import time
import time_local_beta as t


def command(args):
    print('STAGE',args,flush=True)
    subprocess.run(args,check=True,cwd=str(t.w.rt.ROOT))


def main():
    if not (t.OUT/'plan.json').exists():
        records=json.loads((t.OUT/'preparation.json').read_text())['records']
        done=set(int(p.parent.name) for p in t.OUT.glob('parts/*/finalization.json'))
        with mp.Pool(8) as pool:
            while len(done)<len(records):
                ready=[r['source'] for r in records if r['source'] not in done and
                       (r['missing']==0 or (t.OUT/'parts'/('%05d'%r['source'])/'.native_complete').exists())]
                for r in pool.imap_unordered(t.finalize_part,ready):done.add(r['source'])
                print('FINALIZED',len(done),'/',len(records),flush=True)
                if len(done)<len(records):time.sleep(30)
        (t.OUT/'executed_pipeline_source.py').write_bytes(Path(t.__file__).read_bytes())
        command([sys.executable,'analysis/time_local_beta.py','plan','--workers','8'])
    if (t.OUT/'fit_submission.json').exists():
        cluster=json.loads((t.OUT/'fit_submission.json').read_text())['cluster']
    else:
        submitted=subprocess.run(['condor_submit',str(t.OUT/'condor/fit.sub')],check=True,
                    stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
        print(submitted.stdout,flush=True)
        cluster=int(re.search(r'cluster (\d+)',submitted.stdout).group(1))
        t.w.ev.dump(t.OUT/'fit_submission.json',dict(cluster=cluster,complete=False))
    plan=json.loads((t.OUT/'plan.json').read_text())
    folds=plan['folds']
    while True:
        n=sum((t.OUT/'folds'/('%04d'%f['fold'])/'freeze.json').exists() for f in folds)
        print('FROZEN',n,'/',len(folds),flush=True)
        if n==len(folds):break
        query=subprocess.run(['condor_q',str(cluster),'-json'],check=True,stdout=subprocess.PIPE,universal_newlines=True)
        jobs=json.loads(query.stdout.strip() or '[]')
        if not jobs:
            # Queue completion can occur between the filesystem count and query.
            n=sum((t.OUT/'folds'/('%04d'%f['fold'])/'freeze.json').exists() for f in folds)
            if n==len(folds):break
            raise RuntimeError('Training jobs ended with incomplete folds; inspect failure.json')
        held=[j for j in jobs if j['JobStatus']==5]
        if held:raise RuntimeError('Training held: '+str([(j['ProcId'],j.get('HoldReason')) for j in held]))
        time.sleep(30)
    t.w.ev.dump(t.OUT/'fit_submission.json',dict(cluster=cluster,complete=True))
    command([sys.executable,'analysis/time_local_beta.py','evaluate'])
    command([sys.executable,'analysis/verify_time_local_beta.py'])
    print('LOCAL_TIME_EXPERIMENT_COMPLETE',flush=True)


if __name__=='__main__':main()
