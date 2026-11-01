#!/usr/bin/env python3
import subprocess,json,time
from concurrent.futures import ThreadPoolExecutor
from research_joint_adc import OUT,ROOT,CONFIGS
from research_flight_conditional import write
from research_conditional_mc import stamp
jobs=[(dataset,name) for dataset in ['mc','flight'] for name,_,_ in CONFIGS[dataset]]
logs=OUT/'logs';logs.mkdir(exist_ok=True)
def run(row):
 dataset,name=row;path=logs/('%s_%s.log'%row)
 if path.exists():raise FileExistsError(path)
 start=time.monotonic()
 with path.open('w') as stream:
  r=subprocess.run(['python3','analysis/research_joint_adc.py','fit','--dataset',dataset,'--name',name],cwd=str(ROOT),stdout=stream,stderr=subprocess.STDOUT)
 result=dict(dataset=dataset,family=name,exit_code=r.returncode,elapsed_seconds=time.monotonic()-start)
 print('JOINT_ADC_TASK',json.dumps(result),flush=True)
 return result
with ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(run,jobs))
write(OUT/'fit_tasks.json',dict(created_utc=stamp(),complete=all(r['exit_code']==0 for r in rows),tasks=rows))
if any(r['exit_code'] for r in rows):raise RuntimeError('Joint fit failed; preserve and inspect logs')
