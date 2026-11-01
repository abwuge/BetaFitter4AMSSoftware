import json,hashlib
from pathlib import Path
H=Path(__file__).resolve().parent
load=lambda n:json.loads((H/n).read_text())
a=load('native_counter_analysis_summary.json');b=load('extended_summary.json');m=load('mixed_analysis_summary.json')
out={'cohorts':a['cohorts'],'time_convention':{k:{'w68':q['w68'],'gain':q['gain_vs_hybrid_pct']} for k,q in a['results']['R2p5_6'].items() if 'oracle_transport' in k},'lowbeta_covariance':a['covariance']['beta06_09'],'extended_counts':{'events':b['events'],'counts':b['counts']},'extended_models':{k:{'selected_trees':q['selected_trees'],'features':q['features'],'tuning_w68':q['tuning_weighted_w68'],'train_w68':q['training_w68'],'mean_C12':q['test']['C12']['mean'],'mean_C13':q['test']['C13']['mean'],'test_balanced':q['test']['balanced'],'bootstrap':b['bootstrap']['balanced']['methods'][k]} for k,q in b['models'].items()},'extended_mass':b['mass'],'mass_negative':{k:q['balanced_R2p5_6']['negative_mass2_fraction'] for k,q in m['mass'].items()},'mass_covariance':m['mass_linear_covariance']['balanced_R2p5_6'],'table_sha256':hashlib.sha256((H/'report_tables.tsv').read_bytes()).hexdigest()}
(H/'final_check.json').write_text(json.dumps(out,separators=(',',':')))
print('FINAL_CHECK_BYTES',(H/'final_check.json').stat().st_size)
