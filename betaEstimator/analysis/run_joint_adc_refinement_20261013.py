#!/usr/bin/env python3
import subprocess,json,time
from concurrent.futures import ThreadPoolExecutor
from refine_joint_adc_optimizer import OUT,ROOT,REFINE,CONFIGS
from research_flight_conditional import write
from research_conditional_mc import stamp
jobs=[(dataset,name) for dataset in ['mc','flight'] for name,_,_ in CONFIGS[dataset]]
logs=REFINE/'logs';logs.mkdir(exist_ok=True)
def run(row):
 dataset,name=row;path=logs/('%s_%s.log'%row)
 if path.exists():raise FileExistsError(path)
 start=time.monotonic()
 with path.open('w') as stream:
  answer=subprocess.run(['python3','analysis/refine_joint_adc_optimizer.py','fit','--dataset',dataset,'--name',name],cwd=str(ROOT),stdout=stream,stderr=subprocess.STDOUT)
 result=dict(dataset=dataset,family=name,exit_code=answer.returncode,elapsed_seconds=time.monotonic()-start)
 print('JOINT_REFINEMENT_TASK',json.dumps(result),flush=True)
 return result
with ThreadPoolExecutor(max_workers=5) as pool:rows=list(pool.map(run,jobs))
write(REFINE/'fit_tasks.json',dict(complete=all(r['exit_code']==0 for r in rows),created_utc=stamp(),tasks=rows))
if any(r['exit_code'] for r in rows):raise RuntimeError('Joint refinement did not complete; all logs retained')
