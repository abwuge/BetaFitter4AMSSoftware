#!/usr/bin/env python3
"""Compare period-response models using held-out measurements, not beta labels."""
import argparse
import json
from pathlib import Path
import numpy as np
import periodic_joint_adc as periodic
import joint_adc_likelihood as joint
from research_periodic_adc import ROOT, OUT, CONFIGS, context
from research_flight_conditional import load, write, digest
from research_conditional_mc import stamp

DEST = OUT / 'measurement_comparison'


def plan():
    DEST.mkdir(exist_ok=True)
    if (DEST / 'protocol.json').exists():
        raise FileExistsError('Measurement comparison was already declared')
    write(DEST / 'protocol.json', dict(created_utc=stamp(), source_sha256=digest(Path(__file__)),
        fit_protocol_sha256={dataset: digest(OUT / dataset / 'protocol.json') for dataset in ['flight','mc']},
        primary_selection='For each fixed time and rigidity-response family, smallest mean held-out 512-node joint time/ADC negative log likelihood among static, step and spline.',
        comparisons=['step minus static', 'spline minus static'], negative_delta_is_better=True,
        bootstrap=dict(replicates=4096, seed=101320268, unit='entire held-out run, paired across all variants',
            quantiles=[2.5,50,97.5], two_alternative_quantiles=[1.25,50,98.75]),
        uncertainty='Frozen-model evaluation-run variation, not training recalibration or detector-model systematic uncertainty',
        additional_diagnostics='Training-optimum differences, integration convergence, era ADC mean contrasts, fitted noise scales and sensitivity of time/ADC signal response',
        reference_scope='This script reads no RICH or MC truth. Later beta evaluation on the already opened pools is development, not new blind confirmation.',
        production_acceptance=False))
    print('PERIOD_MEASUREMENT_COMPARISON_PROTOCOL', flush=True)


def compare():
    p = json.loads((DEST / 'protocol.json').read_text())
    if digest(Path(__file__)) != p['source_sha256']:
        raise ValueError('Comparison implementation changed')
    if (DEST / 'evaluation.json').exists():
        raise FileExistsError('Comparison already saved')
    result = dict(complete=True, created_utc=stamp(), families=[], RICH_values_read=False, MC_truth_values_read=False,
        protocol_sha256=digest(DEST / 'protocol.json'), selected={})
    rng = np.random.RandomState(p['bootstrap']['seed'])
    for dataset in ['flight','mc']:
        if digest(OUT / dataset / 'protocol.json') != p['fit_protocol_sha256'][dataset]:
            raise ValueError('Scientific fit protocol changed')
        result['selected'][dataset] = {}
        for family, _, _ in CONFIGS[dataset]:
            predictions = []; rows = []
            for variant in periodic.VARIANTS:
                folder = OUT / dataset / (family + '_' + variant)
                fr = json.loads((folder / 'freeze.json').read_text())
                ev = json.loads((folder / 'measurement_evaluation.json').read_text())
                if not ev['complete'] or digest(folder / 'freeze.json') != ev['freeze_sha256'] or digest(folder / 'predictions_before_reference.npz') != ev['predictions_sha256']:
                    raise ValueError('Incomplete or changed period model output')
                v = load(folder / 'predictions_before_reference.npz'); predictions.append(v)
                theta = np.asarray(fr['model']['theta']); nf = (len(theta)-39)//4
                covariance = np.linalg.inv(joint.precision(theta)[1])
                norm = fr['normalization']
                offsets = np.zeros((5,4))
                if variant != 'static':
                    raw = np.eye(5)[:,1:]
                    basis = np.clip((raw-np.asarray(norm['period_center']))/np.asarray(norm['period_scale']), -3,3)
                    offsets = basis.dot(theta[39:].reshape(4,nf)[:,-4:].T)
                    offsets -= offsets[:1]
                start = [json.loads((folder / ('start_%d.json'%i)).read_text()) for i in range(2)]
                rows.append(dict(variant=variant, heldout_nll=float(v['nll'].mean()), material_kappa=theta[:2],
                    lambda_ADC=np.exp(theta[35:39]), period_ADC_log_response_relative_first_era=offsets,
                    period_ADC_gain_relative_first_era=np.exp(offsets), time_scale_diagonal=np.sqrt(np.diag(covariance)[:3]),
                    ADC_scale_diagonal=np.sqrt(np.diag(covariance)[3:]), covariance_eigenvalues=np.linalg.eigvalsh(covariance),
                    training_penalized_nll=fr['model']['penalized_nll'], projected_gradient_max=fr['model']['projected_gradient_max'],
                    start_convergence=[s['numerically_converged'] for s in start],
                    start_training_nll_range=float(np.ptp([s['penalized_nll'] for s in start])),
                    integration_mean256_minus512=ev['quadrature_mean256_minus512'],
                    integration_rms256_minus512=ev['quadrature_rms256_minus512'],
                    time_only_signal_gain_percentiles=ev['time_only_signal_gain_percentiles'],
                    freeze_sha256=digest(folder / 'freeze.json')))
            for v in predictions[1:]:
                for key in ['run','event']:
                    if not np.array_equal(v[key], predictions[0][key]):
                        raise ValueError('Variants were not evaluated on identical events')
            runs = predictions[0]['run']; unique, inverse = np.unique(runs, return_inverse=True)
            nll = np.column_stack([v['nll'] for v in predictions]); delta = nll[:,1:] - nll[:,:1]
            sums = np.column_stack([np.bincount(inverse, weights=delta[:,k]) for k in range(2)])
            counts = np.bincount(inverse); boot = np.empty((p['bootstrap']['replicates'],2))
            for rep in range(len(boot)):
                index = rng.randint(0,len(unique),len(unique)); boot[rep] = sums[index].sum(0)/counts[index].sum()
            comparisons=[]
            for k, variant in enumerate(periodic.VARIANTS[1:]):
                comparisons.append(dict(variant=variant, mean_delta_nll=float(delta[:,k].mean()),
                    paired_run_interval=np.percentile(boot[:,k],p['bootstrap']['quantiles']),
                    two_alternative_interval=np.percentile(boot[:,k],p['bootstrap']['two_alternative_quantiles'])))
            chosen = rows[int(np.argmin([row['heldout_nll'] for row in rows]))]['variant']
            result['selected'][dataset][family] = chosen
            summary = dict(dataset=dataset,family=family,events=len(runs),runs=len(unique),selected=chosen,variants=rows,comparisons=comparisons)
            result['families'].append(summary)
            np.savez_compressed(DEST / (dataset+'_'+family+'_bootstrap.npz'), delta_nll=boot)
            print('PERIOD_MEASUREMENT_RESULT', dataset, family, chosen, json.dumps(comparisons,default=lambda a:a.tolist()), flush=True)
    write(DEST / 'evaluation.json',result)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('stage',choices=['plan','compare']); args=parser.parse_args()
    globals()[args.stage]()
