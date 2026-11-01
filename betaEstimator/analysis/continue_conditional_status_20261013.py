#!/usr/bin/env python3
"""Read-only inventory of inherited continuation; save dated evidence separately."""
import json, hashlib, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/flight_calibration_research'
OUT=BASE/'continuation_20261013'
OUT.mkdir(exist_ok=True)
def digest(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as s:
  for b in iter(lambda:s.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def read(path):
 p=BASE/path
 return json.loads(p.read_text()) if p.exists() else None
result={'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'scope':'Inherited evidence read on 2026-10-13; not counted as new fitting','base':{},'features':[],'domain':{},'protected':{}}
for f in sorted(BASE.rglob('*freeze*.json')):
 if 'continuation_20261013' not in str(f):result['protected'][str(f.relative_to(ROOT))]=digest(f)
for sub in ['conditional_mc','mean_loss_likelihood']:
 for name in ['measured.npz','moments.npz','geometry.root']:
  p=BASE/sub/name
  if p.exists():result['protected'][str(p.relative_to(ROOT))]=digest(p)
for sub in ['scoring/reference_audit.json','scoring/evaluation.json','scoring/joint_bootstrap.json']:
 d=read('conditional_mc/'+sub)
 if d is None:print('MISSING',sub);continue
 if 'reference_audit' in sub:
  result['reference']={k:v for k,v in d.items() if k!='native_source_audits'};print('REFERENCE',json.dumps(result['reference']))
 elif 'evaluation' in sub:
  result['base']['events']=d['events'];result['base']['runs']=d['runs'];result['base']['metrics']={}
  for space in ['u','beta']:
   result['base']['metrics'][space]={}
   for label,row in d['metrics'][space].items():
    result['base']['metrics'][space][label]={'metrics':row['metrics'],'gain':row['versus_ols'],'evaluation_interval':d['paired_evaluation_run_bootstrap'][space][label]}
    if label in ['ols_u','constant_u','native_gaussianR_fitted_u','native_gaussianR_fitted_transport_u','native_studentR5_fitted_u','native_studentR5_fitted_transport_u']:
     print('BASE',space,label,json.dumps(result['base']['metrics'][space][label]))
 else:
  result['base']['joint_intervals']=d['intervals'];result['base']['convergence_failures']=d['convergence_failures'];print('JOINT',d['replicates'],'FAILURES',len(d['convergence_failures']))
for path in sorted((BASE/'feature_noise').glob('*/*/measurement_evaluation.json')):
 d=json.loads(path.read_text());fr=json.loads((path.parent/'freeze.json').read_text())
 row={'dataset':path.parts[-3],'family':path.parts[-2],'measurement':d,'fit':{k:fr['model'][k] for k in ['success','projected_gradient_max','penalized_nll','theta']}}
 result['features'].append(row);print('FEATURE',row['dataset'],row['family'],d['delta_nll'],d['paired_run_delta_nll_interval'],row['fit']['projected_gradient_max'])
d=read('feature_noise/mc/beta_evaluation.json')
if d:
 result['feature_beta']=d
 for row in d['families']:
  for space in ['u','beta']:print('FEATURE_BETA',row['family'],space,json.dumps(row['spaces'][space]))
else:print('MISSING feature beta evaluation')
for name in ['quadrature','propagation']:
 d=read('mean_loss_likelihood/domain_audit/'+name+'.json')
 if d:
  result['domain'][name]={k:v for k,v in d.items() if k!='records'}
  if name=='propagation':
   for row in d['summary']:
    print('PROPAGATION',row['q'],row['kappa'],'valid',row['valid'],'stopped',row['stopped'],'error',row['linear_minus_nonlinear_u'])
  else:
   for row in d['families']:print('QUADRATURE',row['family'],json.dumps({k:v for k,v in row.items() if k not in ['events','family']}))
notes=[json.loads(p.read_text()) for p in (BASE/'conditional_mc/bootstrap').glob('[0-9]*.json')]
result['bootstrap']={'completed':sum(x.get('complete',False) for x in notes),'requested':128}
print('BOOTSTRAP',result['bootstrap'],'PROTECTED',len(result['protected']))
p=OUT/'inherited_status.json'
if p.exists():raise FileExistsError(p)
p.write_text(json.dumps(result,indent=2)+'\n')
print('STATUS',str(p),'SHA256',digest(p))
