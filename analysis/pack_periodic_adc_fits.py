#!/usr/bin/env python3
"""Pack unstarted fits into idle worker capacity; retain original fit protocols."""
import argparse
import json
import os
import socket
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from dispatch_periodic_adc import ROOT, BATCH, RUNTIME
from research_flight_conditional import write,digest
from research_conditional_mc import stamp

PACK = BATCH/'packed'


def query(cluster):
    return json.loads(subprocess.check_output(['condor_q',str(cluster),'-json','-attributes',
        'ClusterId,ProcId,JobStatus,RequestCpus,RemoteHost']).decode() or '[]')


def prepare():
    if PACK.exists():raise FileExistsError('Inspect the existing packed dispatch')
    cluster=json.loads((BATCH/'submission.json').read_text())['cluster']
    before=query(cluster)
    evidence=subprocess.check_output(['condor_status','-constraint','(OpSysAndVer=="CentOS7") && (Cpus>=32)',
        '-af','Machine','State','Activity','Cpus']).decode()
    hosts={}
    for line in evidence.splitlines():
        host,state,activity,cpus=line.split()
        hosts[host]=max(hosts.get(host,0),int(cpus) if state=='Unclaimed' and activity=='Idle' else 0)
    idle=[host for host in hosts if hosts[host]>0]
    if not idle:raise ValueError('No idle suitable worker host')
    workers=min(2,len(idle))
    if len([row for row in before if row['JobStatus']==1])<workers:raise ValueError('Too few idle tasks to pack')
    PACK.mkdir();(PACK/'claims').mkdir()
    write(PACK/'before.json',dict(created_utc=stamp(),queue=before,slot_evidence=evidence,
        reason='Only two suitable physical hosts were idle at the preceding inventory; disjoint thirty-host eligibility leaves most starts waiting. Pack eight four-thread fits per reserved 32-CPU worker.',
        original_dispatch_protocol_sha256=digest(BATCH/'protocol.json'),completed_starts=sorted(str(p.relative_to(ROOT)) for p in BATCH.parent.glob('*/*/start_*.json'))))
    answer=subprocess.check_output(['condor_hold','-constraint','ClusterId==%d && JobStatus==1'%cluster,
        '-reason','Unstarted periodic-ADC fit will run in a recorded packed-capacity worker']).decode()
    after=query(cluster);tasks=sorted(row['ProcId'] for row in after if row['JobStatus']==5 and not row.get('RemoteHost'))
    if not tasks:raise ValueError('No unstarted held tasks are available')
    spec=json.loads((BATCH/'protocol.json').read_text())
    for index in tasks:
        dataset,family,variant,start=spec['tasks'][index]
        path=BATCH.parent/dataset/(family+'_'+variant)/('start_%d.json'%start)
        if path.exists() or (BATCH/('task_%02d.json'%index)).exists():raise ValueError('A held task has completed output; inspect before packing')
    ordered=sorted(hosts,key=lambda host:(-hosts[host],host));pools=[ordered[i::workers] for i in range(workers)]
    worker=PACK/'worker.sh'
    worker.write_text('''#!/usr/bin/env bash
set -eo pipefail
cd {root}
unset PYTHONHOME
export PATH=/usr/bin:/bin:$PATH
export PYTHONPATH={runtime}/python:/home/ams/hxwu/.local/lib
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
/usr/bin/python3 -c 'import numpy,scipy,awkward; assert numpy.__version__=="1.19.1" and scipy.__version__=="0.18.1" and awkward.__version__=="1.9.0"'
exec /usr/bin/python3 analysis/pack_periodic_adc_fits.py worker --index "$1"
'''.format(root=ROOT,runtime=RUNTIME));worker.chmod(0o755)
    submit='\n'.join(['universe = vanilla','executable = '+str(worker),'arguments = $(number)','getenv = true',
        'request_cpus = 32','request_memory = 32GB','should_transfer_files = NO',
        'requirements = (TARGET.OpSysAndVer=="CentOS7") && stringListMember(TARGET.Machine,"$(hosts)",";")',
        'on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)',
        'output = '+str(PACK/'worker_$(number).out'),'error = '+str(PACK/'worker_$(number).err'),'log = '+str(PACK/'condor.log'),
        'queue number,hosts from (']+[str(i)+' '+';'.join(pools[i]) for i in range(workers)]+[')',''])
    (PACK/'fits.sub').write_text(submit)
    write(PACK/'protocol.json',dict(created_utc=stamp(),held_original_cluster=cluster,tasks=tasks,workers=workers,
        subprocesses_per_worker=8,fit_threads=4,reserved_CPUs_per_worker=32,host_pools=pools,
        original_queue_after_hold=after,hold_reply=answer,scientific_protocol_changed=False,
        original_dispatch_source_sha256=digest(ROOT/'analysis/dispatch_periodic_adc.py'),
        source_sha256=digest(Path(__file__)),worker_sha256=digest(worker),submit_sha256=digest(PACK/'fits.sub')))
    reply=subprocess.check_output(['condor_submit','-terse',str(PACK/'fits.sub')]).decode().strip()
    write(PACK/'submission.json',dict(created_utc=stamp(),result=reply,cluster=int(reply.split('.')[0]),protocol_sha256=digest(PACK/'protocol.json')))
    print('PACKED_PERIOD_FITS',len(tasks),'tasks',workers,'workers',reply,flush=True)


def worker(number):
    p=json.loads((PACK/'protocol.json').read_text())
    if digest(Path(__file__))!=p['source_sha256'] or digest(ROOT/'analysis/dispatch_periodic_adc.py')!=p['original_dispatch_source_sha256']:
        raise ValueError('Packed worker source changed')
    def loop(slot):
        done=[]
        for index in p['tasks']:
            claim=PACK/'claims'/('%02d'%index)
            try:claim.mkdir()
            except FileExistsError:continue
            write(claim/'owner.json',dict(created_utc=stamp(),worker=number,slot=slot,host=socket.getfqdn(),pid=os.getpid()))
            with (claim/'fit.log').open('w') as stream:
                r=subprocess.run([sys.executable,'analysis/dispatch_periodic_adc.py','worker','--index',str(index)],cwd=str(ROOT),stdout=stream,stderr=subprocess.STDOUT)
            write(claim/'complete.json',dict(created_utc=stamp(),index=index,exit_code=r.returncode,worker=number,host=socket.getfqdn()))
            if r.returncode:raise RuntimeError('Packed fit %d failed; retain its output'%index)
            done.append(index);print('PACKED_PERIOD_TASK',number,index,flush=True)
        return done
    with ThreadPoolExecutor(max_workers=p['subprocesses_per_worker']) as pool:
        groups=list(pool.map(loop,range(p['subprocesses_per_worker'])))
    write(PACK/('worker_%d.json'%number),dict(complete=True,created_utc=stamp(),host=socket.getfqdn(),tasks=[i for group in groups for i in group]))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','worker']);parser.add_argument('--index',type=int,default=0)
    args=parser.parse_args()
    if args.stage=='worker':worker(args.index)
    else:prepare()
