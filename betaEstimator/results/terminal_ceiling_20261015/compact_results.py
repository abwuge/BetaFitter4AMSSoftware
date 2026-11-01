import json
from pathlib import Path
H=Path(__file__).resolve().parent
b=json.loads((H/'bootstrap_native.json').read_text())
for st in ['all','R2p5_6','beta06_09']:
 print('NATIVE_CI',st,b[st]['n'])
 for name,m in b[st]['methods'].items():
  if name in ['measured_refit','oracle_transport','transport_plus_true_residual_deposit','transport_plus_true_residual']:
   print(name,'gain',round(m['gain_pct'],3),'CI',[round(x,3) for x in m['gain95']],'extra',round(m['additional_gain_pct'],3),[round(x,3) for x in m['additional_gain95']])
a=json.loads((H/'mixed_analysis_summary.json').read_text());print('MIXED_N',a['events'],a['counts'])
for name,r in a['models'].items():
 print('MIXED_MODEL',name,[(st,round(r['test_metrics'][st]['w68'],8),round(r['test_metrics'][st]['weighted_gain_vs_OLS_pct'],3),round(r['test_metrics'][st]['mean'],7)) for st in ['balanced_R2p5_6','C12_R2p5_6','C13_R2p5_6']],r.get('injection'))
for name,r in a['mass'].items():
 print('MASS',name,[(st,round(r[st]['weighted_w68'],6),round(r[st]['positive_mass_relative']['w68'],6)) for st in ['balanced_R2p5_6','C12_R2p5_6','C13_R2p5_6']])
