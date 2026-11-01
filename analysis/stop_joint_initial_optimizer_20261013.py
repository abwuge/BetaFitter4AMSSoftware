#!/usr/bin/env python3
"""Stop only the exact current research runner after retaining completed starts."""
import subprocess,os,signal,time,json
from pathlib import Path
rows=[]
for line in subprocess.check_output(['ps','-u','hxwu','-o','pid=','-o','ppid=','-o','args='],universal_newlines=True).splitlines():
    fields=line.strip().split(None,2)
    if len(fields)==3:rows.append((int(fields[0]),int(fields[1]),fields[2]))
parents=[r for r in rows if r[2]=='python3 analysis/run_joint_adc_fits_20261013.py']
if len(parents)!=1:raise RuntimeError('Expected one exact current research runner: '+repr(parents))
parent=parents[0][0];selected={parent};changed=True
while changed:
    new={pid for pid,ppid,args in rows if ppid in selected};changed=not new.issubset(selected);selected|=new
path=Path('results/flight_calibration_research/continuation_20261013/joint_adc/initial_optimizer_interruption.json')
if path.exists():raise FileExistsError(path)
record=dict(time=time.strftime('%Y-%m-%dT%H:%M:%S%z'),reason='Memory20 unscaled optimizer exhausted 2700-iteration sequences without convergence; numerical preconditioning/memory amendment only, no model or data change',
    parent=parent,processes=[r for r in rows if r[0] in selected],control_tool='Original job control ID no longer existed after a tool-session restart; verified exact current server process lineage before stopping')
path.write_text(json.dumps(record,indent=2)+'\n')
os.kill(parent,signal.SIGTERM)
for pid in sorted(selected-{parent},reverse=True):
    try:os.kill(pid,signal.SIGTERM)
    except ProcessLookupError:pass
print('STOPPED_OWN_JOINT_INITIAL',record,flush=True)
