import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,json,csv
from audit_common import HERE,write_json
from bootstrap_audit import bootstrap
out={};rows=[]
for family,file,names,baseline in [
 ('native','native_counter_analysis_errors.npz',['hybrid','measured_refit','oracle_transport','transport_plus_true_residual_deposit'],'hybrid'),
 ('mixed','mixed_analysis_errors.npz',['ols','measured_tof','beta_fed_tracker','unrestricted','true_rigidity'],'ols'),
 ('extended','extended_errors.npz',['ols','measured_tof','unrestricted_tracker','all_noRich','all_with_RICH','all_noRich_capacity'],'ols')]:
 if not (HERE/file).exists():continue
 a=np.load(str(HERE/file));test=a['test'] if family=='native' else a['role']==1
 scopes=[('all',test),('lowR',test&(a['R']>=2.5)&(a['R']<6))] if family=='native' else [('lowR',test)]
 for st,take in scopes:
  errors={k:a['u_true'][take]/(a['u_true'][take]+a[k][take])-1 for k in names}
  if any(not np.isfinite(v).all() for v in errors.values()):raise ValueError('Nonfinite relative beta')
  b=bootstrap(errors,a['source'][take],baseline,weights=a['weights'][take] if 'weights' in a else None,stratification=a['label'][take] if 'label' in a else None)
  out[family+'_'+st]=b
  for k,m in b['methods'].items():rows.append([family,st,k,b['n'],m['w68'],m['gain_pct'],*m['gain95']])
write_json(HERE/'beta_relative_summary.json',out)
with (HERE/'beta_relative_metrics.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['family','scope','method','n','W68_beta_relative','gain_pct','lower95','upper95']);w.writerows([[format(v,'.10g') if isinstance(v,float) else v for v in r] for r in rows])
print('BETA_RELATIVE_COMPLETE',len(rows),flush=True)
