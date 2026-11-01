#!/usr/bin/env python3
"""Mechanical/provenance verification and scientific-scope summary for round 2."""
import json
import sys
import subprocess
from pathlib import Path
import numpy as np
import research_flight_timing_redundancy as r
ROOT=r.ROOT;OUT=ROOT/'results/flight_calibration_research'
SCRIPTS=['flight_readout_supplement.C','build_flight_readout.sh','research_flight_readout_audit.py',
         'flight_timing_redundancy.py','research_flight_timing_redundancy.py','research_flight_timing_lineage.py',
         'research_newcalib_timing_redundancy.py','research_valid_timing_pipeline.py','research_timing_mc_closure.py',
         'verify_flight_readout_continuation.py']


def get(path):return json.loads(path.read_text())

def main():
    test=subprocess.run([sys.executable,str(ROOT/'tests/test_flight_timing_redundancy.py')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    (OUT/'continuation_unit_tests.log').write_text(test.stdout)
    if test.returncode:raise RuntimeError(test.stdout)
    legacy=get(OUT/'verification.json');assert legacy['passed'] and legacy['old_frozen_models_verified']==999
    audit=get(OUT/'readout_audit/verification.json');protocol=get(OUT/'readout_audit/protocol.json')
    assert audit['complete'] and audit['events']==576 and audit['runs']==80 and audit['association_mismatches']==0
    assert audit['raw_strip_sum_check_passed']
    for row in audit['rows']:
        for field,c in row['comparison'].items():assert c['equal']==c['total'] and c['max_abs']==0,(row['source'],field)
    assert protocol['native_source_sha256']==r.digest(ROOT/'analysis/flight_readout_supplement.C')
    assert protocol['binary_sha256']==r.digest(OUT/'readout_audit/supplement')
    native_summary=dict(events=audit['events'],runs=audit['runs'],association_mismatches=0,
        TOF_comparison='Exact float equality for ADC, LT/HT/SHT raw times, default deposited energy and bar identity on all audited records',
        TOF_default_dynode_hits=sum(x['tof_default_equals_dynode'] for x in audit['rows']),TOF_default_anode_hits=sum(x['tof_default_equals_anode'] for x in audit['rows']),
        TOF_hits=sum(x['tof_positive_hits'] for x in audit['rows']),raw_tracker_strips=sum(x['nstrip'] for x in audit['rows']),
        raw_tracker_overflow_strips=sum(x['overflow_strips'] for x in audit['rows']),nonzero_tracker_strip_status=sum(x['bad_status_strips'] for x in audit['rows']),
        positive_strip_sum_max_abs_difference=max(x['raw_positive_sum_max_abs_difference'] for x in audit['rows']))
    branch_audit=get(OUT/'timing_lineage/audit.json');comparison=get(OUT/'valid_timing_native/baseline_comparison.json')
    time_results={};identities={};calibration_sd={}
    for baseline in ['native','newcalib']:
        folder=OUT/('valid_timing_'+baseline);p=get(folder/'protocol.json');f=get(folder/'freeze.json');e=get(folder/'evaluation.json');b=get(folder/'bootstrap.json')
        assert e['complete'] and b['complete'] and b['requested_replicates']==32 and b['successful_replicates']==32
        assert e['freeze_sha256']==r.digest(folder/'freeze.json')
        assert f['module_sha256']==r.digest(ROOT/'analysis/flight_timing_redundancy.py')
        assert not set(f['response_runs'])&set(f['correction_runs'])
        assert not set(f['training_runs'])&set(f['evaluation_runs'])
        with np.load(str(folder/'measured.npz'),allow_pickle=False) as data:
            identities[baseline]=(data['run'].copy(),data['event'].copy())
            assert np.isfinite(data['t']).all() and not np.any(data['t']==-1) and not np.any(data['t']==0)
            assert len(data['run'])==p['common_events']==142280
            assert not any('rich' in key.lower() or 'truth' in key.lower() for key in data.files)
            cohorts={name:dict(events=int(np.isin(data['fold'],folds).sum()),runs=len(np.unique(data['run'][np.isin(data['fold'],folds)]))) for name,folds in [('response',[1,2]),('correction',[3,4]),('evaluation',[0])]}
        info=dict(cohorts=cohorts,events=e['events'],runs=e['runs'],candidates={})
        for name,result in e['candidates'].items():
            assert result['maximum_signal_response_error']<1.e-10
            info['candidates'][name]=dict(layer_W68_fractional_gain=[1-x['corrected']['w68']/x['baseline']['w68'] for x in result['timing_contrasts']],
                layer_W68_ns=[dict(baseline=x['baseline']['w68'],corrected=x['corrected']['w68']) for x in result['timing_contrasts']],
                layer_RMSE_ns=[dict(baseline=x['baseline']['rmse'],corrected=x['corrected']['rmse']) for x in result['timing_contrasts']],
                inverse_speed_change=result['inverse_speed_change'],signal_response_max_error=result['maximum_signal_response_error'],
                baseline_beta_above_one=result['baseline_beta_above_one'],corrected_beta_above_one=result['corrected_beta_above_one'],
                bootstrap=b['uncertainty'][name],minimum_scaled_projected_Gram_eigenvalue=min(f['candidates'][name]['eigenvalues']))
        time_results[baseline]=info
    for i in range(2):assert np.array_equal(identities['native'][i],identities['newcalib'][i])
    mc=get(OUT/'timing_mc_closure/evaluation.json');mf=get(OUT/'timing_mc_closure/freeze.json')
    assert mc['complete'] and mc['freeze_sha256']==r.digest(OUT/'timing_mc_closure/freeze.json')
    assert not set(mf['train_runs'])&set(mf['test_runs'])
    with np.load(str(OUT/'timing_mc_closure/measured.npz'),allow_pickle=False) as data:
        assert not any('truth' in key.lower() or 'rich' in key.lower() for key in data.files)
        assert not np.any(data['t']==-1)
        mc_counts=dict(events=len(data['run']),runs=len(np.unique(data['run'])))
    mc['measured_training_and_test_counts']=mc_counts;mc['paired_held_run_bootstrap']=get(OUT/'timing_mc_closure/evaluation_bootstrap.json')
    with np.load(str(OUT/'timing_mc_closure/predictions.npz'),allow_pickle=False) as data:
        error=data['baseline_u']-data['truth_effective_u'];decomposition={}
        for name in ['end_difference','all_thresholds']:
            correction=data['baseline_u']-data[name+'_u']
            predicted_change=np.mean(correction*correction)-2*np.mean(error*correction)
            measured_change=np.mean((error-correction)**2)-np.mean(error*error)
            assert abs(predicted_change-measured_change)<1.e-15
            decomposition[name]=dict(correction_mean_square=np.mean(correction**2),twice_error_correction_product=2*np.mean(error*correction),change_in_mean_squared_error=measured_change)
        mc['held_truth_MSE_decomposition']=decomposition
    files=[ROOT/'analysis'/s for s in SCRIPTS]+[ROOT/'tests/test_flight_timing_redundancy.py']
    metadata=dict(status='research_continuation_checkpoint_not_final_algorithm',mechanical_checks_passed=True,new_unit_tests=6,
        old_unit_tests=legacy['unit_test_count'],old_frozen_models_verified=999,readout_audit=native_summary,
        timing_lineage=dict(audited_events=branch_audit['events'],audited_runs=branch_audit['runs'],method0_method1_equal=branch_audit['between_method_equality'],
            original_native_missing_time_events=154,newcalib_missing_or_numerically_invalid_events=106,valid_common_events=142280,valid_common_comparison=comparison),
        flight_timing=time_results,AMS_MC_timing_closure=mc,
        limitations=['Flight scores are projected-time residuals, not beta truth resolution',
            'The response and correction families are low dimensional and do not implement energy-loss transport',
            'Projected Gram rank is not full physical-likelihood information or proof of output identifiability',
            'Run-held samples have been used for development diagnostics; no new RICH sample was opened',
            'MC truth target is effective TOF timing slope, not L5 beta or RICH beta',
            'No MC coefficients or truth-fitted parameters were transferred to flight',
            '32 full-run flight refits describe sampling variation, not all detector-model uncertainties',
            'Signal injection preserves supplied geometry/ADC and valid raw-time status, but does not rerun the complete upstream beta-fed reconstruction',
            'Unexplained timing/geometry tails remain in the no-residual-cut flight sample; RMSE must not be suppressed',
            'Absolute energy scale, passive energy losses, physical mixture/rigidity likelihood and new blind RICH beta evaluation are still open'],
        preserved_failures=['Initial native input list lacked reused-event runs; fixed by union of both existing native lists',
            'First coherent injection changed missing zero thresholds into apparent signals; fixed injection, same frozen coefficients',
            'V1 finite-only layer-time validity admitted -1; retain V1 diagnostics and use only validity-corrected V2 for main comparison',
            'V1 NewCalib bootstrap cancelled to avoid continued work on invalid input; its partial progress is not a completed interval'],
        script_sha256={str(p.relative_to(ROOT)):r.digest(p) for p in files})
    r.write(OUT/'continuation_summary.json',metadata)
    compact=dict(passed=True,new_tests=6,old_tests=legacy['unit_test_count'],old_models_unchanged=999,native_event_matches=576,
        common_valid_flight_events=142280,AMS_MC_holdout_events=mc['events'],RICH_evaluation=False,
        scope='Mechanical provenance, numerical identities and defined experiments; scientific improvement not certified',summary_sha256=r.digest(OUT/'continuation_summary.json'))
    r.write(OUT/'continuation_verification.json',compact)
    print(json.dumps(r.native(compact),indent=2),flush=True)

if __name__=='__main__':main()
