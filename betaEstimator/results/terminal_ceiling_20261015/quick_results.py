import json
from pathlib import Path
H=Path(__file__).resolve().parent
p=H/'native_counter_analysis_summary.json'
a=json.loads(p.read_text())
print('NATIVE_COHORTS',a['cohorts'])
for st in ['all_beta_gt06','R2p5_6','beta06_09','beta09_1']:
 print('NATIVE_STRATUM',st)
 for name,m in a['results'][st].items():
  print(name,m['n'],round(m['w68'],8),round(m['mean'],8),round(m['std'],8),round(m['gain_vs_hybrid_pct'],3))
 print('COV',a['covariance'][st])
 if st in a['geometry']:print('GEOM',{k:v for k,v in a['geometry'][st].items() if k!='covariance'})
p=H/'mixed_analysis_summary.json'
if p.exists():
 a=json.loads(p.read_text());print('MIXED_N',a['events'],a['counts'])
 for name,r in a['models'].items():
  print('MIXED_MODEL',name,[(st,round(r['test_metrics'][st]['w68'],8),round(r['test_metrics'][st]['weighted_gain_vs_OLS_pct'],3)) for st in ['balanced_R2p5_6','C12_R2p5_6','C13_R2p5_6']],r.get('injection'))
 for name,r in a['mass'].items():
  print('MASS',name,[(st,r[st]['weighted_w68'],r[st]['positive_mass_relative']['w68']) for st in ['balanced_R2p5_6','C12_R2p5_6','C13_R2p5_6']])
