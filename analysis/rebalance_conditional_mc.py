#!/usr/bin/env python3
"""Return completed workers' disjoint host groups to queued MC workers."""
import argparse
import json
import subprocess
import time
from pathlib import Path
from research_conditional_mc import OUT,stamp
from research_flight_conditional import write

BATCH=OUT/'batch'
CLUSTER=397520


def rebalance():
    state_path=BATCH/'host_ownership.json'
    if state_path.exists():state=json.loads(state_path.read_text())
    else:
        submission=json.loads((BATCH/'submission.json').read_text())
        state=dict(cluster=CLUSTER,host_owner={host:i for i,hosts in enumerate(submission['disjoint_groups']) for host in hosts},history=[])
    if state['cluster']!=CLUSTER:raise ValueError('Unexpected queue cluster')
    result=subprocess.check_output(['condor_q',str(CLUSTER),'-json','-attributes','ProcId,JobStatus,RemoteHost,NumJobStarts'],universal_newlines=True)
    jobs=json.loads(result);active={r['ProcId'] for r in jobs};idle=[r['ProcId'] for r in jobs if r['JobStatus']==1 and not r.get('RemoteHost')]
    running=[r for r in jobs if r['JobStatus']==2];running_hosts={r['RemoteHost'].split('@')[-1] for r in running if r.get('RemoteHost')}
    if len(running_hosts)!=len(running):raise ValueError('Workers share a physical host')
    donor=[host for host,owner in state['host_owner'].items() if owner not in active and host not in running_hosts]
    changes=[]
    if idle and donor:
        # One recipient owns all released hosts until it completes. No existing
        # running worker or its candidate group is moved.
        recipient=min(idle,key=lambda i:(sum(v==i for v in state['host_owner'].values()),i))
        for host in donor:state['host_owner'][host]=recipient
        hosts=sorted(h for h,o in state['host_owner'].items() if o==recipient)
        expression='(TARGET.OpSysAndVer == "CentOS7") && stringListMember(TARGET.Machine,"'+ ';'.join(hosts)+'",";") && (TARGET.Arch == "X86_64") && (TARGET.Disk >= RequestDisk) && (TARGET.Memory >= RequestMemory) && (TARGET.Cpus >= RequestCpus) && (TARGET.FileSystemDomain == MY.FileSystemDomain)'
        constraint='ClusterId == %d && ProcId == %d && JobStatus == 1 && RemoteHost =?= UNDEFINED'%(CLUSTER,recipient)
        answer=subprocess.check_output(['condor_qedit','-constraint',constraint,'Requirements',expression],universal_newlines=True)
        changes.append(dict(recipient=recipient,released_hosts=donor,result=answer.strip()))
    note=dict(created_utc=stamp(),jobs=len(jobs),running=len(running),distinct_running_hosts=len(running_hosts),idle=len(idle),changes=changes)
    state['history'].append(note);write(state_path,state)
    print('HOST_ALLOCATION',json.dumps(note),flush=True)
    return len(jobs)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--watch',action='store_true');args=parser.parse_args()
    while True:
        count=rebalance()
        if not args.watch or not count:break
        time.sleep(30)
