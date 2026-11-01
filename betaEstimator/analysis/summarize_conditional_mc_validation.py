#!/usr/bin/env python3
"""Create a compact, evidence-linked view of the ongoing continuation."""
import json
import subprocess
from pathlib import Path
from research_conditional_mc import OUT,ROOT,stamp
from research_flight_conditional import write,digest


def read(path):return json.loads(path.read_text()) if path.exists() else None


def main():
    result=dict(created_utc=stamp(),scope='Continuation progress; no global optimality or impossibility claim',sources={},base_models=[],feature_models=[],domain=[])
    refpath=OUT/'scoring/reference_audit.json';ref=read(refpath)
    if ref:
        result['reference']={k:v for k,v in ref.items() if k not in ['native_source_audits']};result['sources'][str(refpath.relative_to(ROOT))]=digest(refpath)
    path=OUT/'scoring/evaluation.json';base=read(path)
    if base:
        for label,record in base['metrics']['u'].items():
            beta=base['metrics']['beta'][label]
            row=dict(label=label,u_metrics=record['metrics'],u_gain=record['versus_ols'],beta_metrics=beta['metrics'],beta_gain=beta['versus_ols'],
                u_interval=base['paired_evaluation_run_bootstrap']['u'][label],beta_interval=base['paired_evaluation_run_bootstrap']['beta'][label])
            result['base_models'].append(row)
        result['sources'][str(path.relative_to(ROOT))]=digest(path)
    path=ROOT/'results/flight_calibration_research/feature_noise/mc/beta_evaluation.json';feature=read(path)
    if feature:
        result['feature_beta']=feature;result['sources'][str(path.relative_to(ROOT))]=digest(path)
    for path in (ROOT/'results/flight_calibration_research/feature_noise').glob('*/*/measurement_evaluation.json'):
        row=read(path);row.update(dataset=path.parts[-3],family=path.parts[-2]);result['feature_models'].append(row);result['sources'][str(path.relative_to(ROOT))]=digest(path)
    path=ROOT/'results/flight_calibration_research/mean_loss_likelihood/domain_audit/propagation.json';domain=read(path)
    if domain:
        for row in domain['summary']:
            error=row['linear_minus_nonlinear_u'];result['domain'].append(dict(q=row['q'],kappa=row['kappa'],tested=row['configurations'],valid=row['valid'],stopped=row['stopped'],
                u_error_rms=error['rms'] if error else None,u_error_range=[error['percentiles'][0],error['percentiles'][-1]] if error else None))
        result['sources'][str(path.relative_to(ROOT))]=digest(path)
    files=sorted((OUT/'bootstrap').glob('[0-9]*.json'));notes=[read(path) for path in files]
    result['bootstrap_progress']=dict(requested=read(OUT/'bootstrap/protocol.json')['replicates'],completed=sum(bool(r.get('complete')) for r in notes),
        failed=[path.name for path,r in zip(files,notes) if not r.get('complete')],
        elapsed_seconds=[r['elapsed_seconds'] for r in notes if 'elapsed_seconds' in r])
    joint=read(OUT/'scoring/joint_bootstrap.json')
    if joint:result['base_joint_intervals']=joint['intervals'];result['base_joint_convergence_failures']=joint['convergence_failures']
    result['open_items']=['New independent flight RICH resolution validation has not been opened','C13 primary raw reference remains outside the audited set','Detector response and event-selection assumptions are not proven by numerical checks']
    write(OUT/'validation_progress.json',result)
    if ref:print('REFERENCE',ref['accepted'],ref['rejection_counts'])
    for row in result['base_models']:
        if row['label'] in ['ols_u','native_gaussianR_fitted_u','native_gaussianR_fitted_transport_u','native_studentR5_fitted_u','native_studentR5_fitted_transport_u']:
            print('BASE',row['label'],'u gain',row['u_gain']['w68_improvement'],'interval',row['u_interval']['w68_improvement'],'beta gain',row['beta_gain']['w68_improvement'],'interval',row['beta_interval']['w68_improvement'])
    if feature:
        for row in feature['families']:
            for space in ['u','beta']:
                cc=row['spaces'][space]['comparisons']['ols'];print('FEATURE',row['family'],space,'W68 gain',cc['point']['w68_improvement'],'CI',cc['interval']['w68_improvement'],'RMSE gain',cc['point']['rmse_improvement'])
    for row in result['feature_models']:print('FEATURE_LIKELIHOOD',row['dataset'],row['family'],row['delta_nll'],row['paired_run_delta_nll_interval'])
    print('BOOTSTRAP',result['bootstrap_progress']['completed'],'of',result['bootstrap_progress']['requested'])
    for row in result['domain']:
        if row['q']>=.4:print('LATENT_DOMAIN',row)

if __name__=='__main__':main()
