#!/usr/bin/env python3
"""Compact, text-only report evidence. Never exports raw events or models."""
import json,hashlib,csv,subprocess
from pathlib import Path
H=Path(__file__).resolve().parent
load=lambda n:json.loads((H/n).read_text())
def keep(x,fields):return {k:x[k] for k in fields if k in x}
def met(x):return keep(x,['n','mean','median','w68','std','rmse','tail02','tail05'])
def compact_boot(x):return {k:keep(v,['w68','gain_pct','gain95','additional_gain_pct','additional_gain95']) for k,v in x['methods'].items()}
E={'scope':'MC audit performed in situ on ssh:sdiat; compact scalar summaries only. No large-file transfer.','directory':str(H),'reference_date':'2026-10-15','metrics':'u=1/beta; W68=(q84-q16)/2; gain=100*(1-W68_new/W68_reference). All confidence intervals are paired whole-source bootstrap conditional on fitted models.','learner_test':load('learner_self_test.json')}
E['broad_counts']={str(z):keep(load('broad_Z%d_manifest.json'%z),['events','train','test','previous_comparison']) for z in [1,2,6,8]}
E['broad']={}
for z in [1,2,6,8]:
 rows=[];boot=load('bootstrap_broad_Z%d.json'%z)
 for p in sorted(H.glob('fit_Z%d_*.json'%z)):
  a=load(p.name);row=keep(a,['tag','group','feature_count','train_available','train','tune','test_events','paired_test','selected_trees','loss'])
  for key,st in [('all','all_beta_gt06'),('lowR','R2p5_6')]:
   m=a['test'][st];bb=boot['all' if key=='all' else 'R2p5_6']['methods'][a['tag']]
   row[key]={'model':met(m),'hybrid':met(m['previous_hybridBeta']),'gain_pct':m['gain_vs_hybrid_pct'],'gain95':bb['gain95']}
  row['injection']=a['injection'];row['tuning_width']=[v['w68'] for v in a['tuning_metrics']]
  rows.append(row)
 E['broad'][str(z)]=rows
E['native_collection']=load('native_collection_summary.json')
E['counter_provenance']=load('counter_plane_provenance.json')
E['native_readout_verification']=load('native_readout_verification.json')
a=load('native_counter_analysis_summary.json');boot=load('bootstrap_native.json')
E['native']={'cohorts':a['cohorts'],'strata':{}}
for st,bkey in [('all_beta_gt06','all'),('R2p5_6','R2p5_6'),('beta06_09','beta06_09')]:
 r={}
 for name,m in a['results'][st].items():
  r[name]=dict(met(m),gain_pct=m['gain_vs_hybrid_pct'])
  if name in boot[bkey]['methods']:r[name].update(keep(boot[bkey]['methods'][name],['gain95','additional_gain_pct','additional_gain95']))
 E['native']['strata'][st]={'methods':r,'covariance':a['covariance'][st],'geometry':a['geometry'].get(st)}
E['native_lowR_focus']=load('native_lowR_focus_summary.json')
a=load('mixed_analysis_summary.json');boot=load('bootstrap_mixed.json')
E['mixed']={'events':a['events'],'counts':a['counts'],'models':{},'mass':{},'prior_tests':a['prior_weight_tests'],'caveats':[a['reference_caveat'],a['domain_caveat']]}
for name,m in a['models'].items():
 row=keep(m,['features','train','tune','test','selected_trees','injection'])
 row['strata']={}
 for st,bkey in [('balanced_R2p5_6','balanced'),('C12_R2p5_6','C12'),('C13_R2p5_6','C13')]:
  q=m['test_metrics'][st];val=dict(met(q),weighted_w68=q['weighted_w68'],gain_pct=q['weighted_gain_vs_OLS_pct'])
  if name in boot[bkey]['methods']:val.update(keep(boot[bkey]['methods'][name],['gain95','additional_gain_pct','additional_gain95']))
  row['strata'][bkey]=val
 E['mixed']['models'][name]=row
massboot=load('bootstrap_mass.json')
for name,m in a['mass'].items():
 row={}
 for st,bkey in [('balanced_R2p5_6','balanced'),('C12_R2p5_6','C12'),('C13_R2p5_6','C13')]:
  q=m[st];val=dict(met(q),weighted_w68=q['weighted_w68'],negative_fraction=q['negative_mass2_fraction'],positive_mass_w68=q['positive_mass_relative']['w68'])
  if name in massboot[bkey]['methods']:val.update(keep(massboot[bkey]['methods'][name],['gain_pct','gain95']))
  row[bkey]=val
 E['mixed']['mass'][name]=row
E['mixed_mass_covariance']=a['mass_linear_covariance']['balanced_R2p5_6']
E['mixed_manifest']=keep(load('mixed_manifest.json'),['events','primary_valid','mass_caveat'])
a=load('extended_summary.json');E['extended']={'events':a['events'],'counts':a['counts'],'claim_scope':a['claim_scope'],'models':{},'mass':a['mass'],'plan':load('extended_plan.json')}
for name,m in a['models'].items():
 row=keep(m,['features','selected_trees','tuning_weighted_w68','training_w68']);row['strata']={}
 for st,q in m['test'].items():row['strata'][st]=dict(met(q),weighted_w68=q['weighted_w68'],gain_pct=q['gain_pct'],**keep(a['bootstrap'][st]['methods'][name],['gain95','additional_gain_pct','additional_gain95']))
 E['extended']['models'][name]=row
E['gaussian']=load('gaussian_surrogate_summary.json')
E['beta_bin_profiles']={}
for z in [1,2,6,8]:
 a=load('fit_Z%d_measured_all_t640_l15_m160_f100_core.json'%z)
 E['beta_bin_profiles'][str(z)]={k:{'n':v['n'],'new_w68':v['w68'],'hybrid_w68':v['previous_hybridBeta']['w68'],'gain_pct':v['gain_vs_hybrid_pct']} for k,v in a['test'].items() if k.startswith('beta_')}
files=sorted([p for p in H.iterdir() if p.is_file() and p.suffix in ['.py','.cpp','.C','.sh','.json'] and p.name!='report_evidence.json'])
E['source_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files if p.suffix!='.json'}
with (H/'audit_text_manifest.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['file','bytes','sha256']);w.writerows((p.name,p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest()) for p in files)
(H/'report_evidence.json').write_text(json.dumps(E,ensure_ascii=False,separators=(',',':'),allow_nan=False))
print('EVIDENCE_COMPLETE',len((H/'report_evidence.json').read_bytes()),'files',len(files),flush=True)
