#!/usr/bin/env python3
"""Verify the layer-use repair and write the final continuation result index."""
import json
import sys
import subprocess
import numpy as np
import research_flight_timing_redundancy as r
ROOT=r.ROOT;BASE=ROOT/'results/flight_calibration_research'


def read(p):return json.loads(p.read_text())

def main():
    audit=read(BASE/'used_layer_audit/flight_audit.json');prior=read(BASE/'continuation_summary.json');flight={};identities=[]
    for name in ['native','newcalib']:
        folder=BASE/('four_layer_timing_'+name);protocol=read(folder/'protocol.json');freeze=read(folder/'freeze.json');ev=read(folder/'evaluation.json');boot=read(folder/'bootstrap.json')
        assert ev['complete'] and boot['complete'] and boot['successful_replicates']==boot['requested_replicates']==32
        assert ev['freeze_sha256']==r.digest(folder/'freeze.json')
        assert freeze['module_sha256']==r.digest(ROOT/'analysis/flight_timing_redundancy.py')
        assert not set(freeze['response_runs'])&set(freeze['correction_runs'])
        assert not set(freeze['training_runs'])&set(freeze['evaluation_runs'])
        with np.load(str(folder/'measured.npz'),allow_pickle=False) as data:
            assert not any('rich' in k.lower() or 'truth' in k.lower() for k in data.files)
            identities.append((data['run'].copy(),data['event'].copy()))
            assert len(data['run'])==audit['kept_four_used_events']==protocol['events']
            assert np.isfinite(data['t']).all() and not np.any(data['t']==-1) and not np.any(data['t']==0)
            cohorts={label:dict(events=int(np.isin(data['fold'],folds).sum()),runs=len(np.unique(data['run'][np.isin(data['fold'],folds)]))) for label,folds in [('response',[1,2]),('correction',[3,4]),('evaluation',[0])]}
        info=dict(cohorts=cohorts,events=ev['events'],runs=ev['runs'],candidates={})
        for candidate,v in ev['candidates'].items():
            assert v['maximum_signal_response_error']<1.e-10
            info['candidates'][candidate]=dict(layer_widths_ns=[dict(baseline=x['baseline']['w68'],corrected=x['corrected']['w68']) for x in v['timing_contrasts']],
                layer_W68_gain_percent=[100*(1-x['corrected']['w68']/x['baseline']['w68']) for x in v['timing_contrasts']],
                layer_RMSE_ns=[dict(baseline=x['baseline']['rmse'],corrected=x['corrected']['rmse']) for x in v['timing_contrasts']],
                signal_response_error=v['maximum_signal_response_error'],u_change=v['inverse_speed_change'],
                beta_above_one_before=v['baseline_beta_above_one'],beta_above_one_after=v['corrected_beta_above_one'],
                refit_bootstrap=boot['uncertainty'][candidate])
        flight[name]=info
    for i in range(2):assert np.array_equal(identities[0][i],identities[1][i])
    with np.load(str(BASE/'used_layer_audit/flight_quality.npz'),allow_pickle=False) as quality:
        wanted=(quality['tof_hsumhu']==4)&(quality['tof_hsumh']==4)
        assert np.array_equal(identities[0][0],quality['run'][wanted]) and np.array_equal(identities[0][1],quality['event'][wanted])
    folder=BASE/'four_layer_timing_mc_closure';mc=read(folder/'evaluation.json');mf=read(folder/'freeze.json');mp=read(folder/'protocol.json');mb=read(folder/'evaluation_bootstrap.json')
    assert mc['freeze_sha256']==r.digest(folder/'freeze.json') and not set(mf['train_runs'])&set(mf['test_runs'])
    assert mf['module_sha256']==r.digest(ROOT/'analysis/flight_timing_redundancy.py')
    with np.load(str(folder/'measured.npz'),allow_pickle=False) as a:
        assert not any('truth' in k.lower() or 'rich' in k.lower() for k in a.files)
        mc['prepared_events']=len(a['run']);mc['prepared_runs']=len(np.unique(a['run']))
    mc['four_used_selection']=mp;mc['paired_held_run_bootstrap']=mb
    with np.load(str(folder/'predictions.npz'),allow_pickle=False) as a:
        e=a['baseline_u']-a['truth_effective_u'];mc['mse_decomposition']={}
        for candidate in ['end_difference','all_thresholds']:
            c=a['baseline_u']-a[candidate+'_u'];change=float(np.mean(c*c)-2*np.mean(e*c))
            assert abs(change-(np.mean((e-c)**2)-np.mean(e**2)))<1.e-15
            mc['mse_decomposition'][candidate]=dict(correction_mean_square=float(np.mean(c*c)),twice_error_correction_product=float(2*np.mean(e*c)),change_in_MSE=change)
    verified=0
    for p in (ROOT/'results/time_local_beta/folds').glob('*/freeze.json'):
        for block,h in read(p)['model_sha256'].items():
            assert r.digest(p.parent/('block_%s_model.json'%block))==h;verified+=1
    assert verified==999
    for test in ['test_flight_calibration_core.py','test_flight_timing_redundancy.py']:
        result=subprocess.run([sys.executable,str(ROOT/'tests'/test)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
        (BASE/('final_'+test+'.log')).write_text(result.stdout)
        assert result.returncode==0,result.stdout
    out=dict(status='context_boundary_research_handoff_not_final_algorithm',readout_audit=prior['readout_audit'],
        timing_lineage=prior['timing_lineage'],used_layer_audit=audit,flight_timing=flight,AMS_MC_timing_closure=mc,
        prior_result_index='continuation_summary.json records V1/V2 studies before used-layer repair and must not be substituted for this index',
        current_result_index='four_layer_continuation_summary.json',
        tests=dict(new=6,old=6,old_frozen_models_unchanged=verified),
        limitations=['No deployed flight energy-loss model, no MeV absolute calibration, no passive-loss reconstruction',
            'Flight score is projected timing contrast, not beta resolution',
            'Only native four-used-layer domain tested; native hit selection is not independently unbiased noise sampling',
            'No RICH measurement opened; no new blind validation claimed',
            'MC target is stored effective TOF truth-time slope, not primary L5 beta or RICH beta',
            'MC parameters are isolated and never imported into flight',
            '32 flight whole-run refits assess sampling, not full detector model uncertainty; 512 MC resamples keep estimator frozen',
            'Only the supplied raw-readout/time interface was coherently injected; upstream beta-fed reconstruction not rerun',
            'Unchanged simplified response/correction family is not the complete physical energy-loss likelihood',
            'No proof of global optimum or impossibility of the broader project'],
        result_sha256={str(p.relative_to(BASE)):r.digest(p) for p in [BASE/'readout_audit/verification.json',BASE/'used_layer_audit/flight_audit.json']+
            [BASE/name/file for name in ['four_layer_timing_native','four_layer_timing_newcalib','four_layer_timing_mc_closure'] for file in ['protocol.json','freeze.json','evaluation.json']]},
        script_sha256=dict(prior['script_sha256'],**{str(p.relative_to(ROOT)):r.digest(p) for p in [ROOT/'analysis/research_four_layer_timing.py',ROOT/'analysis/verify_four_layer_timing.py']}))
    r.write(BASE/'four_layer_continuation_summary.json',out)
    record=dict(passed=True,scope='Provenance, used-layer domain, tests and unchanged protected models only; not scientific acceptance',
        events=audit['kept_four_used_events'],flight_evaluation_events=flight['native']['events'],MC_evaluation_events=mc['events'],
        old_models_unchanged=verified,summary_sha256=r.digest(BASE/'four_layer_continuation_summary.json'))
    r.write(BASE/'four_layer_verification.json',record);print(json.dumps(record,indent=2),flush=True)

if __name__=='__main__':main()
