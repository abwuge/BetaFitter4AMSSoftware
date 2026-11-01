#!/usr/bin/env python3
"""Reassign only completed workers' host groups to idle workers of our cluster."""
import json,subprocess
from pathlib import Path
from research_amplitude_noise_20261013 import OUT
from research_flight_conditional import write
from research_conditional_mc import stamp
BATCH=OUT/'bootstrap/batch';CLUSTER=398913
submission=json.loads((BATCH/'submission.json').read_text())
if submission['cluster']!=CLUSTER:raise ValueError('Unexpected current research cluster')
statepath=BATCH/'rebalance_state.json'
if statepath.exists():state=json.loads(statepath.read_text())
else:state=dict(cluster=CLUSTER,owner={h:i for i,g in enumerate(submission['disjoint_groups']) for h in g},history=[])
if state['cluster']!=CLUSTER:raise ValueError('Wrong ownership state')
jobs=json.loads(subprocess.check_output(['condor_q',str(CLUSTER),'-json','-attributes','ProcId,JobStatus,RemoteHost'],universal_newlines=True))
active={r['ProcId'] for r in jobs};idle=[r['ProcId'] for r in jobs if r['JobStatus']==1 and not r.get('RemoteHost')]
running=[r for r in jobs if r['JobStatus']==2];hosts={r['RemoteHost'].split('@')[-1] for r in running}
if len(hosts)!=len(running):raise ValueError('Physical-host separation violated')
donors=sorted(h for h,owner in state['owner'].items() if owner not in active and h not in hosts)
changes=[]
if idle and donors:
    for h in donors:
        recipient=min(idle,key=lambda i:(sum(v==i for v in state['owner'].values()),i))
        state['owner'][h]=recipient
    for i in idle:
        owned=sorted(h for h,v in state['owner'].items() if v==i)
        expression='(TARGET.OpSysAndVer == "CentOS7") && stringListMember(TARGET.Machine,"'+ ';'.join(owned)+'",";") && (TARGET.Arch == "X86_64") && (TARGET.Disk >= RequestDisk) && (TARGET.Memory >= RequestMemory) && (TARGET.Cpus >= RequestCpus) && (TARGET.FileSystemDomain == MY.FileSystemDomain)'
        constraint='ClusterId == %d && ProcId == %d && JobStatus == 1 && RemoteHost =?= UNDEFINED'%(CLUSTER,i)
        answer=subprocess.check_output(['condor_qedit','-constraint',constraint,'Requirements',expression],universal_newlines=True)
        changes.append(dict(worker=i,hosts=owned,response=answer.strip()))
record=dict(created_utc=stamp(),jobs=len(jobs),running=len(running),idle=len(idle),released_hosts=donors,changes=changes)
state['history'].append(record);write(statepath,state)
print('AMPLITUDE_REBALANCE',json.dumps(dict(jobs=len(jobs),running=len(running),idle=len(idle),released=len(donors),updated=len(changes))),flush=True)
