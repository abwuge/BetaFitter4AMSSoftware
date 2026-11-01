import json,subprocess,time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
HERE=Path(__file__).resolve().parent
plan=[]
for z in [1,2,8]:
 for group in ['geometry','geometry_time','measured_tof','measured_all']:
  plan.append(dict(z=z,group=group,trees=640 if group=='measured_all' else 160,leaves=15,minleaf=160,fraction=1.,loss='core'))
for group in ['raw_readout','truth_deposit','beta_fed','unrestricted','rich_fusion','shuffle']:
 plan.append(dict(z=6,group=group,trees=160 if group=='shuffle' else 640,leaves=15,minleaf=160,fraction=1.,loss='core'))
plan.append(dict(z=6,group='measured_all',trees=1280,leaves=31,minleaf=80,fraction=1.,loss='core'))
plan.append(dict(z=6,group='measured_all',trees=640,leaves=15,minleaf=160,fraction=1.,loss='mse'))
for fraction in [.25,.5]:plan.append(dict(z=6,group='measured_all',trees=640,leaves=15,minleaf=160,fraction=fraction,loss='core'))
(HERE/'fit_queue_plan.json').write_text(json.dumps(dict(created=time.time(),selection='Only tuning sources select checkpoint; whole testing sources untouched by parameter fitting.',plan=plan),indent=2))
def run(x):
 tag='Z%d_%s_t%d_l%d_m%d_f%d_%s'%(x['z'],x['group'],x['trees'],x['leaves'],x['minleaf'],round(100*x['fraction']),x['loss'])
 dest=HERE/('fit_'+tag+'.json')
 if dest.exists():return dict(tag=tag,status='already_complete')
 cmd=['python3',str(HERE/'train_broad.py')]
 for k,v in x.items():cmd+=['--'+k,str(v)]
 start=time.monotonic()
 with (HERE/('queue_'+tag+'.log')).open('w') as log:r=subprocess.run(cmd,cwd=str(HERE),stdout=log,stderr=subprocess.STDOUT)
 ans=dict(tag=tag,exit_code=r.returncode,complete=dest.exists(),seconds=time.monotonic()-start)
 print('FIT_QUEUE',ans,flush=True)
 return ans
with ThreadPoolExecutor(max_workers=8) as pool:results=list(pool.map(run,plan))
(HERE/'fit_queue_complete.json').write_text(json.dumps(results,indent=2))
if any(not x.get('complete',x.get('status')=='already_complete') for x in results):raise RuntimeError('One or more fits failed')
