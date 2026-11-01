#!/usr/bin/env python3
"""Collect this continuation's evidence without rewriting older result indexes."""
import datetime
import hashlib
import json
import subprocess
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/flight_calibration_research'
OUT=BASE/'mean_loss_likelihood'

def read(p):return json.loads(Path(p).read_text())
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    audit=read(BASE/'native_timing_audit/verification.json');geom=read(OUT/'geometry_verification.json');prop=read(OUT/'propagation_verification.json');fit=read(OUT/'contrast_diagnostic.json')
    assert audit['complete'] and audit['events']==98 and audit['association_mismatches']==0
    for row in audit['rows']:
        assert all(row['identity'].values())
        for key in ['stored_time_vs_native_cluster','stored_length_vs_native_beta_length','native_beta_time_vs_cluster_time','stored_chi2_vs_native','tof_rtime','tof_rtime_ht','tof_rtime_sht']:
            assert row['comparisons'][key]['max_abs']==0
        assert all(row['extreme']['used']) and all(row['extreme']['goodtime'])
    assert geom['events']==geom['valid']==12000 and geom['max_interlayer_path_over_c_error_ns']<2e-6
    assert digest(OUT/'geometry.root')==geom['geometry_sha256'] and digest(OUT/'moments.npz')==geom['moments_sha256']
    for name,h in prop['source_sha256'].items():assert digest(ROOT/'analysis'/name)==h
    assert prop['complete'] and len(prop['cases'])==8 and prop['zero_loss_time_ns']['max_abs']<1e-10
    assert fit['complete'] and fit['status']=='diagnostic_only_not_a_deployable_calibration'
    assert len(fit['models'])==len(fit['evaluation'])==len(fit['bootstrap'])==4
    for row in fit['evaluation']:assert digest(OUT/(row['tag']+'_contrast_freeze.json'))==row['freeze_sha256']
    assert all(len(row['coefficient_samples'])==64 for row in fit['bootstrap'])
    with np.load(OUT/'measured.npz',allow_pickle=False) as f:
        assert not any('rich' in k.lower() or 'truth' in k.lower() or k.startswith('mc_') for k in f.files)
        ids=(f['run'].astype(np.uint64)<<np.uint64(32))|f['event'].astype(np.uint64)
        assert len(np.unique(ids))==12000
        assert not np.intersect1d(f['run'][f['fold']!=0],f['run'][f['fold']==0]).size
    conditional=OUT/'conditional';numeric=read(conditional/'numerical_verification.json');outputs=read(conditional/'output_diagnostic.json')
    assert numeric['passed'] and outputs['complete'] and len(outputs['records'])==3
    for name,h in numeric['source_sha256'].items():assert digest(ROOT/'analysis'/name)==h
    assert digest(conditional/'libconditional.so')==numeric['library_sha256']
    assert digest(ROOT/'analysis/research_conditional_output.py')==outputs['script_sha256']
    assert abs(numeric['density_integrals']['48']-1)<2e-5
    conditional_fits=[]
    for row in outputs['records']:
        folder=conditional/row['label'];protocol=read(folder/'protocol.json');frozen=read(folder/'freeze.json');ev=read(folder/'evaluation.json')
        assert len(frozen['starts'])==3 and ev['complete']
        assert digest(folder/'freeze.json')==ev['freeze_sha256']==row['freeze_sha256']
        assert digest(ROOT/'analysis/flight_conditional_likelihood.py')==frozen['module_sha256']
        assert digest(ROOT/'analysis/flight_conditional_likelihood.cpp')==frozen['cpp_sha256']
        assert digest(conditional/'libconditional.so')==frozen['library_sha256']
        assert protocol['events']==9356 and ev['events']==2644 and not protocol['RICH_read'] and not protocol['MC_read']
        for variant in row['free_slope_output'].values():assert variant['max_injected_response_error']<1e-10
        conditional_fits.append(dict(label=row['label'],protocol=protocol,model=frozen['model'],
            start_nlls=[s['nll'] for s in frozen['starts']],zero_loss_nll=frozen['null_model']['nll'],unit_material_nll=frozen['unit_material_model']['nll'],
            evaluation=ev,output_diagnostic=row))
    verified=0
    for p in (ROOT/'results/time_local_beta/folds').glob('*/freeze.json'):
        for block,h in read(p)['model_sha256'].items():
            assert digest(p.parent/('block_%s_model.json'%block))==h;verified+=1
    assert verified==999
    plan=read(ROOT/'results/time_local_beta/plan.json')
    assert digest(ROOT/'results/time_local_beta/training.npz')==plan['training_sha256']
    tests={}
    for name in ['test_flight_calibration_core.py','test_flight_timing_redundancy.py','test_flight_mean_loss.py']:
        r=subprocess.run([sys.executable,str(ROOT/'tests'/name)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
        (OUT/(name+'.log')).write_text(r.stdout)
        if r.returncode:raise RuntimeError(r.stdout)
        tests[name]=dict(passed=True,log_sha256=digest(OUT/(name+'.log')))
    sources=['flight_native_timing_audit.C','verify_flight_native_timing_audit.py','flight_material_geometry.C','build_flight_material.sh',
             'research_flight_material.py','flight_mean_loss.py','verify_flight_mean_loss.py','research_mean_loss_contrasts.py','verify_mean_loss_continuation.py',
             'build_flight_native_timing_audit.sh','flight_conditional_likelihood.cpp','flight_conditional_likelihood.py',
             'build_flight_conditional.sh','research_flight_conditional.py','research_conditional_output.py']
    artifact_paths=[BASE/'native_timing_audit/protocol.json',BASE/'native_timing_audit/verification.json',BASE/'native_timing_audit/training_membership.json',
                    OUT/'protocol.json',OUT/'geometry_verification.json',OUT/'propagation_verification.json',OUT/'contrast_diagnostic.json',
                    conditional/'numerical_verification.json',conditional/'output_diagnostic.json']
    for row in outputs['records']:
        artifact_paths.extend(conditional/row['label']/name for name in ['protocol.json','freeze.json','evaluation.json'])
    summary=dict(created_at=datetime.datetime.now().astimezone().isoformat(),status='research_continuation_not_final_algorithm',
        prior_index='four_layer_continuation_summary.json; prior experiments are not new work of this continuation',
        native_timing_audit=dict(events=98,association_mismatches=0,all_stored_times_lengths_raw_threshold_times_match_exactly=True,
            extremes=[dict(run=r['extreme']['run'],event=r['extreme']['event'],native_beta=r['extreme']['beta_value'],native_chi2=r['extreme']['chi2'],
                           four_used=True,four_good_time=True,pulse_counts=r['extreme']['pulse_count']) for r in audit['rows']],
            final_training_membership=audit['training_membership'],
            interpretation='Not ntuple transcription or chosen track/BetaH mismatch. Multiple readout pulses are available; causal pulse assignment remains unresolved.'),
        native_material_geometry=geom,
        propagation=dict(tested_events=128,cases=prop['cases'],zero_loss_time_ns=prop['zero_loss_time_ns'],small_loss_derivative_time_ns=prop['small_loss_derivative_time_ns'],
            scope=prop['scope'],physical_forward_model_only=True,measured_beta_above_one_inference_not_implemented=True),
        contrast_diagnostic={k:fit[k] for k in ['status','train_events','train_runs','evaluation_events','evaluation_runs','models','evaluation','bootstrap','omissions','box_note']},
        conditional_working_likelihood=dict(numerical_checks=numeric,fits=conditional_fits,
            native_noise_family_free_slope_output_difference=outputs['native_rigidity_noise_family_weight_output_difference'],
            boundary_status='First implemented normalized conditional timing model with true-rigidity integration and fitted isotope/spectrum nuisance; not the complete selected-readout likelihood',
            output_status='Free-slope timing diagnostic implemented without event isotope/spectrum shrinkage; not yet independently validated for beta resolution'),
        checks=dict(unit_tests=tests,unit_test_count=19,conditional_gradient_dimensions=13,conditional_rigidity_noise_variants_checked=2,
            conditional_normalization_checked=True,old_frozen_models_unchanged=verified,old_final_training_cache_unchanged=True,RICH_read=False,MC_read=False),
        not_completed=['Full normalized selected raw-readout likelihood with actual native hit selection and ADC response','Native timing-calibration covariance and charge/epoch conditioning',
            'Conditional-model whole-run calibration refits and latent-domain/nonlinear-transport sensitivity','Absolute visible deposit response and passive straggling calibration',
            'New independent AMS MC mean-loss refit and L5/RICH closure','New blind flight RICH evaluation','Global optimality or impossibility proof'],
        source_sha256={str(Path('analysis')/name):digest(ROOT/'analysis'/name) for name in sources},
        result_sha256={str(p.relative_to(BASE)):digest(p) for p in artifact_paths})
    (BASE/'material_continuation_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    record=dict(passed=True,scope='Numerical/provenance tests and frozen-file integrity; the physical calibration model is not accepted',
        unit_test_count=19,conditional_fits_verified=3,conditional_starts_verified=9,conditional_density_and_gradients_verified=True,
        old_models_unchanged=verified,old_final_training_cache_unchanged=True,
        summary_sha256=digest(BASE/'material_continuation_summary.json'),source_sha256=digest(__file__))
    (BASE/'material_continuation_verification.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2),flush=True)

if __name__=='__main__':main()
