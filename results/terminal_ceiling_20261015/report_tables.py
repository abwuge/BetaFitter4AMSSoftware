import json,csv,hashlib
from pathlib import Path
H=Path(__file__).resolve().parent
load=lambda n:json.loads((H/n).read_text())
rows=[]
def add(group,scope,name,n,m,w=None,gain=None,ci=None,extra=None):
 if ci is None:ci=m.get('gain95',['',''])
 extra=extra or m
 rows.append([group,scope,name,n,w if w is not None else m.get('w68',''),gain if gain is not None else m.get('gain_pct',''),*ci,m.get('mean',''),m.get('std',''),extra.get('additional_gain_pct',''),*extra.get('additional_gain95',['',''])])
for z in [1,2,6,8]:
 a=load('fit_Z%d_measured_all_t640_l15_m160_f100_core.json'%z);boot=load('bootstrap_broad_Z%d.json'%z)
 for scope,st,bkey in [('all','all_beta_gt06','all'),('lowR','R2p5_6','R2p5_6')]:
  m=a['test'][st];q=boot[bkey]['methods'][a['tag']];add('broad',scope,str(z),m['n'],m,gain=m['gain_vs_hybrid_pct'],ci=q['gain95'])
for p in sorted(H.glob('fit_Z6_*.json')):
 a=load(p.name);b=load('bootstrap_broad_Z6.json')['all']['methods'][a['tag']];m=a['test']['all_beta_gt06'];add('C12_ablation','all',a['tag'],m['n'],m,gain=m['gain_vs_hybrid_pct'],ci=b['gain95'])
a=load('native_counter_analysis_summary.json');b=load('bootstrap_native.json')
for scope,st,bkey in [('all','all_beta_gt06','all'),('lowR','R2p5_6','R2p5_6')]:
 for name in ['hybrid','measured_refit','oracle_transport','transport_plus_measured','transport_plus_true_residual','transport_plus_true_residual_deposit','oracle_noise_removed','all_truth']:
  m=a['results'][st][name];q=b[bkey]['methods'].get(name,{});add('native',scope,name,m['n'],m,gain=m['gain_vs_hybrid_pct'],ci=q.get('gain95',['','']),extra=q)
a=load('native_lowR_focus_summary.json')['bootstrap']
for name in ['selected_measured','selected_transport_residual','selected_transport_residual_deposit']:add('native_focus','lowR',name,a['n'],a['methods'][name])
a=load('mixed_analysis_summary.json');boot=load('bootstrap_mixed.json')
for name,q in a['models'].items():
 m=q['test_metrics']['balanced_R2p5_6'];b=boot['balanced']['methods'].get(name,{});add('mixed','lowR_balanced',name,m['n'],m,w=m['weighted_w68'],gain=m['weighted_gain_vs_OLS_pct'],ci=b.get('gain95',['','']),extra=b)
a=load('mixed_analysis_summary.json');boot=load('bootstrap_mass.json')
for name,q in a['mass'].items():
 m=q['balanced_R2p5_6'];b=boot['balanced']['methods'].get(name,{});add('mass2','lowR_balanced',name,m['n'],m,w=m['weighted_w68'],gain=b.get('gain_pct',0 if name=='OLS_beta_measured_R' else 100 if name=='all_truth' else ''),ci=b.get('gain95',['','']))
if (H/'extended_summary.json').exists():
 a=load('extended_summary.json')
 for name,q in a['models'].items():
  m=q['test']['balanced'];b=a['bootstrap']['balanced']['methods'][name];add('extended','lowR_balanced',name,m['n'],m,w=m['weighted_w68'],gain=m['gain_pct'],ci=b['gain95'],extra=b)
with (H/'report_tables.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['group','scope','method','n','w68','gain_pct','lower95','upper95','mean','std','extra_gain_pct','extra_lower95','extra_upper95'])
 for row in rows:w.writerow([format(x,'.10g') if isinstance(x,float) else x for x in row])
notes={}
notes['counts']={str(z):{k:load('broad_Z%d_manifest.json'%z)[k] for k in ['events','train','test']} for z in [1,2,6,8]}
notes['native_collection']=load('native_collection_summary.json')
notes['code_hashes']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in H.iterdir() if p.suffix in ['.py','.cpp','.C','.sh']}
notes['learner_self_test']=load('learner_self_test.json')
notes['mixed_counts']=load('mixed_analysis_summary.json')['counts']
notes['mixed_species_bias']={name:{st:{k:q[k] for k in ['n','mean','std','w68','weighted_gain_vs_OLS_pct']} for st,q in m['test_metrics'].items() if st in ['C12_R2p5_6','C13_R2p5_6']} for name,m in load('mixed_analysis_summary.json')['models'].items() if name in ['measured_tof','unrestricted','true_rigidity']}
notes['gaussian']={st:{k:v for k,v in q.items() if k not in ['covariances','cell_counts','assumptions']} for st,q in load('gaussian_surrogate_summary.json')['results'].items()}
notes['native_lowR_plan']=load('native_lowR_plan.json')
notes['broad_fit_counts']={str(z):{k:load('fit_Z%d_measured_all_t640_l15_m160_f100_core.json'%z).get(k) for k in ['train','tune','test_events','paired_test','selected_trees']} for z in [1,2,6,8]}
(H/'report_notes.json').write_text(json.dumps(notes,separators=(',',':')))
print('TABLES',len(rows),(H/'report_tables.tsv').stat().st_size,'NOTES',(H/'report_notes.json').stat().st_size)
