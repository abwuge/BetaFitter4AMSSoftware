#!/usr/bin/env python3
"""Submit disjoint physical-host groups; no existing scheduler jobs are edited."""
import json,subprocess,re
from pathlib import Path
from research_amplitude_noise_20261013 import ROOT,OUT
from research_flight_conditional import write,digest
from research_conditional_mc import stamp
BATCH=OUT/'bootstrap/batch'
BATCH.mkdir(exist_ok=True)
if (BATCH/'submission.json').exists():raise FileExistsError('Already submitted')
a=json.loads(subprocess.check_output(['condor_status','-json','-attributes','Machine,OpSysAndVer,State,Activity,Cpus,Memory'],universal_newlines=True))
allhosts=sorted(set(x['Machine'] for x in a if x['OpSysAndVer']=='CentOS7'))
free=sorted(set(x['Machine'] for x in a if x['OpSysAndVer']=='CentOS7' and x['State']=='Unclaimed'))
ordered=free+[h for h in allhosts if h not in free]
workers=32
if len(ordered)<workers:raise ValueError('Insufficient distinct eligible hosts')
groups=[ordered[i::workers] for i in range(workers)]
assert len(set(h for g in groups for h in g))==sum(map(len,groups))
text='''universe = vanilla
executable = %s
arguments = $(worker)
getenv = true
request_cpus = 32
request_memory = 8GB
requirements = (TARGET.OpSysAndVer == "CentOS7") && stringListMember(TARGET.Machine, "$(allowed_hosts)", ";")
should_transfer_files = NO
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
output = %s/worker_$(Process).out
error = %s/worker_$(Process).err
log = %s/condor.log
queue worker,allowed_hosts from (
'''%(ROOT/'analysis/run_amplitude_bootstrap_20261013.sh',BATCH,BATCH,BATCH)
text+='\n'.join('%d %s'%(i,';'.join(g)) for i,g in enumerate(groups))+'\n)\n'
path=BATCH/'bootstrap.sub';path.write_text(text)
answer=subprocess.check_output(['condor_submit',str(path)],universal_newlines=True)
found=re.search(r'cluster\s+(\d+)',answer)
if not found:raise RuntimeError('Submission returned no cluster identifier: '+answer)
record=dict(created_utc=stamp(),cluster=int(found.group(1)),workers=workers,replicates=128,candidate_hosts=ordered,initially_unclaimed_hosts=free,
    disjoint_groups=groups,simultaneous_replicates_per_worker=4,OpenMP_threads_per_replicate=4,
    request_cpus_reason='Existing CentOS7 static-slot START requires RequestCpus == 32; 4 concurrent replica processes each use 4 OpenMP threads',
    submit_sha256=digest(path),wrapper_sha256=digest(ROOT/'analysis/run_amplitude_bootstrap_20261013.sh'),submission_output=answer)
write(BATCH/'submission.json',record)
print('AMPLITUDE_BOOTSTRAP_SUBMITTED',json.dumps(record),flush=True)
