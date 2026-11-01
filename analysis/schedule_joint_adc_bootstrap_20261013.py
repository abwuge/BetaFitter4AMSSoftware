#!/usr/bin/env python3
"""Submit and rebalance only the current joint-ADC bootstrap cluster."""
import argparse,json,subprocess,re
from bootstrap_joint_adc_20261013 import OUT,ROOT
from research_flight_conditional import write,digest
from research_conditional_mc import stamp
BATCH=OUT/'batch'


def submit():
    BATCH.mkdir(exist_ok=True)
    if (BATCH/'submission.json').exists():raise FileExistsError('Joint bootstrap already submitted')
    protocol=json.loads((OUT/'protocol.json').read_text())
    if protocol['replicates']!=128:raise ValueError('Wrapper requires declared 128 replicates')
    queue_text=subprocess.check_output(['condor_q','hxwu','-json','-attributes','ProcId,ClusterId,JobStatus,RemoteHost'],universal_newlines=True)
    own=json.loads(queue_text) if queue_text.strip() else []
    occupied={r['RemoteHost'].split('@')[-1] for r in own if r.get('RemoteHost') and r['JobStatus']==2}
    slots=json.loads(subprocess.check_output(['condor_status','-json','-attributes','Machine,OpSysAndVer,State'],universal_newlines=True))
    hosts=sorted({r['Machine'] for r in slots if r['OpSysAndVer']=='CentOS7'}-occupied)
    free=sorted({r['Machine'] for r in slots if r['OpSysAndVer']=='CentOS7' and r['State']=='Unclaimed'}-occupied)
    ordered=free+[h for h in hosts if h not in free];workers=32
    if len(ordered)<workers:raise ValueError('Insufficient disjoint physical hosts')
    groups=[ordered[i::workers] for i in range(workers)]
    wrapper=BATCH/'worker.sh'
    wrapper.write_text('''#!/usr/bin/env bash
set -eo pipefail
cd %s
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
worker="$1"
pids=()
for offset in 0 32 64 96; do
  rep=$((worker+offset))
  python3 analysis/bootstrap_joint_adc_20261013.py worker --replicate "$rep" > "%s/rep_${rep}.log" 2>&1 &
  pids+=("$!")
done
status=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then status=1; fi
done
printf 'JOINT_BOOTSTRAP_WORKER_DONE worker=%%s status=%%s host=%%s\\n' "$worker" "$status" "$(hostname)"
exit "$status"
'''%(ROOT,BATCH));wrapper.chmod(0o755)
    text='''universe = vanilla
executable = %s
arguments = $(worker)
getenv = true
request_cpus = 32
request_memory = 8GB
requirements = (TARGET.OpSysAndVer == "CentOS7") && stringListMember(TARGET.Machine,"$(hosts)",";")
should_transfer_files = NO
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
output = %s/worker_$(Process).out
error = %s/worker_$(Process).err
log = %s/condor.log
queue worker,hosts from (
'''%(wrapper,BATCH,BATCH,BATCH)
    text+='\n'.join('%d %s'%(i,';'.join(group)) for i,group in enumerate(groups))+'\n)\n'
    path=BATCH/'bootstrap.sub';path.write_text(text);response=subprocess.check_output(['condor_submit',str(path)],universal_newlines=True)
    match=re.search(r'cluster\s+(\d+)',response)
    if not match:raise RuntimeError('No scheduler cluster ID in submission response '+response)
    record=dict(created_utc=stamp(),cluster=int(match.group(1)),workers=workers,replicates=128,configs=len(protocol['records']),
        excluded_preexisting_own_hosts=sorted(occupied),disjoint_groups=groups,actual_OpenMP_threads_per_host=16,
        request_cpus_reason='Existing CentOS7 static-slot START requires RequestCpus32; four independent repetitions each use four OpenMP threads',
        protocol_sha256=digest(OUT/'protocol.json'),wrapper_sha256=digest(wrapper),submission_sha256=digest(path),submission_response=response)
    write(BATCH/'submission.json',record);print('JOINT_BOOTSTRAP_SUBMITTED',record['cluster'],len(ordered),'hosts',flush=True)


def rebalance():
    submission=json.loads((BATCH/'submission.json').read_text());cluster=submission['cluster'];path=BATCH/'rebalance_state.json'
    state=json.loads(path.read_text()) if path.exists() else dict(cluster=cluster,owner={h:i for i,g in enumerate(submission['disjoint_groups']) for h in g},history=[])
    if state['cluster']!=cluster:raise ValueError('Wrong scheduler ownership state')
    queue_text=subprocess.check_output(['condor_q',str(cluster),'-json','-attributes','ProcId,JobStatus,RemoteHost,HoldReason'],universal_newlines=True)
    jobs=json.loads(queue_text) if queue_text.strip() else []
    active={r['ProcId'] for r in jobs};idle=[r['ProcId'] for r in jobs if r['JobStatus']==1 and not r.get('RemoteHost')]
    running=[r for r in jobs if r['JobStatus']==2];hosts={r['RemoteHost'].split('@')[-1] for r in running}
    if len(hosts)!=len(running):raise ValueError('Joint bootstrap physical-host overlap')
    donors=sorted(h for h,owner in state['owner'].items() if owner not in active and h not in hosts);changes=[]
    if idle and donors:
        for h in donors:
            receiver=min(idle,key=lambda i:(sum(v==i for v in state['owner'].values()),i));state['owner'][h]=receiver
        for i in idle:
            owned=sorted(h for h,v in state['owner'].items() if v==i)
            expression='(TARGET.OpSysAndVer == "CentOS7") && stringListMember(TARGET.Machine,"'+';'.join(owned)+'",";") && (TARGET.Arch == "X86_64") && (TARGET.Disk >= RequestDisk) && (TARGET.Memory >= RequestMemory) && (TARGET.Cpus >= RequestCpus) && (TARGET.FileSystemDomain == MY.FileSystemDomain)'
            constraint='ClusterId == %d && ProcId == %d && JobStatus == 1 && RemoteHost =?= UNDEFINED'%(cluster,i)
            response=subprocess.check_output(['condor_qedit','-constraint',constraint,'Requirements',expression],universal_newlines=True)
            changes.append(dict(worker=i,hosts=owned,response=response.strip()))
    state['history'].append(dict(created_utc=stamp(),running=len(running),idle=len(idle),donor_hosts=donors,changes=changes))
    write(path,state);print('JOINT_BOOTSTRAP_STATE',cluster,'running',len(running),'idle',len(idle),'updated',len(changes),'held',[r for r in jobs if r['JobStatus']==5],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['submit','rebalance']);args=p.parse_args();globals()[args.stage]()
