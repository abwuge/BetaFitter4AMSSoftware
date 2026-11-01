#!/usr/bin/env python3
"""Eight independent conditional-noise starts, sharing available worker capacity."""
import argparse
import json
import os
import socket
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from research_conditional_noise import ROOT,OUT
from research_flight_conditional import write,digest
from research_conditional_mc import stamp
from dispatch_periodic_adc import RUNTIME
BATCH=OUT/'batch'


def prepare():
    if BATCH.exists():raise FileExistsError('Conditional-noise dispatch already exists')
    verified=json.loads((OUT/'numerical_verification.json').read_text())
    if not verified['complete'] or verified['protocol_sha256']!=digest(OUT/'protocol.json'):raise ValueError('Verify the noise implementation first')
    queue=json.loads(subprocess.check_output(['condor_q','hxwu','-json','-attributes','ClusterId,ProcId,JobStatus,RemoteHost']).decode() or '[]')
    if queue:raise ValueError('Inspect existing user jobs before dispatch')
    text=subprocess.check_output(['condor_status','-constraint','(OpSysAndVer=="CentOS7") && (Cpus>=32)',
        '-af','Machine','State','Activity','Cpus']).decode()
    hosts={}
    for line in text.splitlines():
        host,state,activity,cpus=line.split();hosts[host]=max(hosts.get(host,0),int(cpus) if state=='Unclaimed' and activity=='Idle' else 0)
    free=[host for host in hosts if hosts[host]>0]
    if not free:raise ValueError('No worker capacity is currently available')
    workers=min(8,len(free));parallel=(8+workers-1)//workers
    ordered=sorted(hosts,key=lambda host:(-hosts[host],host));pools=[ordered[i::workers] for i in range(workers)]
    tasks=[(dataset,variant,start) for dataset in ['flight','mc'] for variant in ['geometry','geometry_period'] for start in range(2)]
    BATCH.mkdir();(BATCH/'claims').mkdir();worker=BATCH/'worker.sh'
    worker.write_text('''#!/usr/bin/env bash
set -eo pipefail
cd {root}
unset PYTHONHOME
export PATH=/usr/bin:/bin:$PATH
export PYTHONPATH={runtime}/python:/home/ams/hxwu/.local/lib
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
/usr/bin/python3 -c 'import numpy,scipy,awkward; assert numpy.__version__=="1.19.1" and scipy.__version__=="0.18.1" and awkward.__version__=="1.9.0"'
exec /usr/bin/python3 analysis/dispatch_conditional_noise.py worker --index "$1"
'''.format(root=ROOT,runtime=RUNTIME));worker.chmod(0o755)
    submission='\n'.join(['universe = vanilla','executable = '+str(worker),'arguments = $(number)','getenv = true',
        'request_cpus = 32','request_memory = 32GB','should_transfer_files = NO',
        'requirements = (TARGET.OpSysAndVer=="CentOS7") && stringListMember(TARGET.Machine,"$(hosts)",";")',
        'on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)',
        'output = '+str(BATCH/'worker_$(number).out'),'error = '+str(BATCH/'worker_$(number).err'),'log = '+str(BATCH/'condor.log'),
        'queue number,hosts from (']+[str(i)+' '+';'.join(pools[i]) for i in range(workers)]+[')',''])
    (BATCH/'fits.sub').write_text(submission)
    write(BATCH/'protocol.json',dict(created_utc=stamp(),tasks=tasks,workers=workers,parallel_fits_per_worker=parallel,
        threads_per_fit=4,physical_host_pools=pools,available_host_inventory=text,
        runtime_manifest_sha256=digest(RUNTIME/'manifest.json'),scientific_protocol_sha256=digest(OUT/'protocol.json'),
        numerical_verification_sha256=digest(OUT/'numerical_verification.json'),source_sha256=digest(Path(__file__)),
        worker_sha256=digest(worker),submit_sha256=digest(BATCH/'fits.sub'),
        policy='At most eight fit tasks; use distinct physical-host pools where available, packing four-thread fits within reserved 32-CPU workers when capacity is smaller.'))
    answer=subprocess.check_output(['condor_submit','-terse',str(BATCH/'fits.sub')]).decode().strip()
    write(BATCH/'submission.json',dict(created_utc=stamp(),result=answer,cluster=int(answer.split('.')[0]),protocol_sha256=digest(BATCH/'protocol.json')))
    print('CONDITIONAL_NOISE_SUBMITTED',len(tasks),'starts',workers,'workers',answer,flush=True)


def worker(number):
    p=json.loads((BATCH/'protocol.json').read_text())
    if digest(Path(__file__))!=p['source_sha256'] or digest(OUT/'protocol.json')!=p['scientific_protocol_sha256']:
        raise ValueError('Noise worker implementation or scientific protocol changed')
    def one(slot):
        completed=[]
        for index,(dataset,variant,start) in enumerate(p['tasks']):
            claim=BATCH/'claims'/('%02d'%index)
            try:claim.mkdir()
            except FileExistsError:continue
            result_path=OUT/dataset/variant/('start_%d.json'%start)
            if result_path.exists():raise FileExistsError('A noise fit result already exists')
            begin=time.monotonic();write(claim/'owner.json',dict(created_utc=stamp(),task=p['tasks'][index],worker=number,slot=slot,host=socket.getfqdn(),pid=os.getpid()))
            with (claim/'fit.log').open('w') as stream:
                result=subprocess.run([sys.executable,'analysis/research_conditional_noise.py','start','--dataset',dataset,'--variant',variant,'--start',str(start)],cwd=str(ROOT),stdout=stream,stderr=subprocess.STDOUT)
            receipt=dict(created_utc=stamp(),task=p['tasks'][index],index=index,host=socket.getfqdn(),exit_code=result.returncode,elapsed_seconds=time.monotonic()-begin)
            if result_path.exists():
                receipt['result_sha256']=digest(result_path);receipt['numerically_converged']=json.loads(result_path.read_text())['numerically_converged']
            write(claim/'complete.json',receipt)
            if result.returncode:raise RuntimeError('Noise fit failed; retain its trace '+str(claim))
            completed.append(index);print('CONDITIONAL_NOISE_TASK',index,receipt.get('numerically_converged'),flush=True)
        return completed
    with ThreadPoolExecutor(max_workers=p['parallel_fits_per_worker']) as pool:groups=list(pool.map(one,range(p['parallel_fits_per_worker'])))
    write(BATCH/('worker_%d.json'%number),dict(complete=True,created_utc=stamp(),host=socket.getfqdn(),completed_tasks=[index for group in groups for index in group]))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','worker']);parser.add_argument('--index',type=int,default=0);args=parser.parse_args()
    if args.stage=='worker':worker(args.index)
    else:prepare()
