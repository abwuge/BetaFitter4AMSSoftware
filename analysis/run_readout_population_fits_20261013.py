#!/usr/bin/env python3
import subprocess,json,time
from concurrent.futures import ThreadPoolExecutor
from research_readout_population import OUT,ROOT
from research_flight_conditional import write
from research_conditional_mc import stamp
jobs=[('mc',name,variant) for name in ['native_gaussianR','native_studentR5'] for variant in ['population_only','joint']]
jobs += [('flight',name,'joint') for name in ['native_gaussianR','newcalib_gaussianR','native_studentR5']]
logs=OUT/'logs';logs.mkdir(exist_ok=True)
def run(row):
 dataset,name,variant=row;path=logs/('%s_%s_%s.log'%row)
 if path.exists():raise FileExistsError(path)
 start=time.monotonic()
 with path.open('w') as stream:
  r=subprocess.run(['python3','analysis/research_readout_population.py','fit','--dataset',dataset,'--name',name,'--variant',variant],cwd=str(ROOT),stdout=stream,stderr=subprocess.STDOUT)
 result=dict(dataset=dataset,family=name,variant=variant,exit_code=r.returncode,elapsed_seconds=time.monotonic()-start)
 print('POPULATION_TASK',json.dumps(result),flush=True)
 return result
with ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(run,jobs))
write(OUT/'fit_tasks.json',dict(created_utc=stamp(),complete=all(r['exit_code']==0 for r in rows),tasks=rows))
if any(r['exit_code'] for r in rows):raise RuntimeError('Population fit failed; inspect retained starts/logs')
