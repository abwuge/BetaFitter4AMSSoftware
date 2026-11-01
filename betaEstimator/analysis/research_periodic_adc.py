#!/usr/bin/env python3
"""Fit declared epoch ADC nuisances using the original measurement-only pools."""
import argparse
import json
from pathlib import Path
import numpy as np
import joint_adc_likelihood as original
import flight_conditional_likelihood as timing
import free_timing_transport as transport
import periodic_joint_adc as periodic
from research_amplitude_noise_20261013 import inputs, CONFIGS
from joint_adc_scaled_optimizer import solve_scaled
from research_flight_conditional import load, write, digest
from research_conditional_mc import stamp

ROOT = periodic.ROOT
OUT = periodic.OUT
SOURCES = ['periodic_joint_adc.py', 'research_periodic_adc.py', 'joint_adc_likelihood.py',
           'joint_adc_likelihood.cpp', 'free_timing_transport.py', 'joint_adc_scaled_optimizer.py']


def prepare(dataset):
    OUT.mkdir(exist_ok=True); destination = OUT / dataset; destination.mkdir(exist_ok=True)
    if (destination / 'protocol.json').exists():
        raise FileExistsError('Period-response protocol already frozen')
    folder, a, g, adc, domain, train, test = inputs(dataset)
    variants = {}
    for variant in periodic.VARIANTS:
        feature, norm = periodic.features(a, train, dataset, variant)
        variants[variant] = dict(normalization=norm, nfeatures=feature.shape[2],
            nparameters=39 + 4 * feature.shape[2], bounds=original.bounds(feature.shape[2]))
    seeds = {}
    for family, field, nu in CONFIGS[dataset]:
        path = original.OUT / dataset / family / 'freeze.json'
        seeds[family] = dict(path=str(path), sha256=digest(path))
    write(destination / 'protocol.json', dict(created_utc=stamp(), dataset=dataset, variants=variants,
        configs=CONFIGS[dataset], seeds=seeds, train_events=int(train.sum()), evaluation_events=int(test.sum()),
        train_runs=np.unique(a['run'][train]), evaluation_runs=np.unique(a['run'][test]),
        hypothesis='Shared ADC noise may be inflated by time-dependent light/gain response. Add only measured-era ADC intercepts; never fit an additive TOF time shift or per-event isotope posterior.',
        ADC_mean='Existing geometry/bar response plus four declared era contrast features per layer, plus exp(log_lambda)*(u_L5-1.15). Layer speed slopes and the noise scale remain common across eras.',
        flight_eras=periodic.YEARS, period_variants=['No epoch dependence', 'Fixed UTC-era offsets', 'Continuous linear interpolation between fixed era midpoint gain anchors'],
        MC_control='The MC epoch coordinate is a deterministic run hash, explicitly a negative control for unsupported run-dependent fitting; it is not evidence about actual flight-era detector changes.',
        nodes=256, evaluation_nodes=512, penalty=.01,
        penalty_formula='0.5 * 0.01 * sum(squared non-intercept ADC response coefficients), including standardized era coefficients',
        starts=2, maximum_rounds=3, projected_gradient_max=1e-3,
        optimizer='Existing measured-curvature-scaled L-BFGS-B, maximum 6000 iterations per round; unchanged model, bounds and objective',
        initializations=['Own previously frozen joint model padded with zero epoch coefficients', 'Existing coefficient-free physical timing initialization, lambda=3 and own ADC medians'],
        start_selection='Minimum penalized training NLL among the two completed starts; numerical convergence must be checked and failures are retained',
        model_selection='Within each fixed time/rigidity-response family, select the smallest held-out measurement joint NLL among the three declared variants. No RICH residual or MC truth enters this choice.',
        uncertainty_and_scope='A new independent flight reference cohort is required for confirmatory performance. The already opened 6751-event prospective sample is now development evidence and must not be called blind again.',
        resolution_interpretation='Separate changes in period-dependent bias from within-period beta width; a narrower pooled distribution due only to bias removal is not proof of lower stochastic beta error.',
        original_training_pool_only=True, RICH_targets_used_in_fit=False, MC_truth_targets_used_in_fit=False,
        previously_opened_RICH_is_development=True, MC_coefficients_transferred_to_flight=False,
        upstream_timing_calibration_reopened=False, old_production_models_replaced=False,
        original_prospective_verification_sha256=digest(ROOT / 'results/flight_calibration_research/continuation_20261013/prospective_flight/verification.json'),
        measurement_sha256=digest(folder / 'measured.npz'), geometry_sha256=digest(folder / 'moments.npz'),
        source_sha256={name: digest(ROOT / 'analysis' / name) for name in SOURCES},
        unchanged_joint_library_sha256=digest(original.OUT / 'libjointadc.so')))
    print('PERIOD_ADC_PROTOCOL', dataset, int(train.sum()), int(test.sum()), len(CONFIGS[dataset]) * len(periodic.VARIANTS), 'models', flush=True)


def context(dataset, family, variant):
    parent = OUT / dataset; p = json.loads((parent / 'protocol.json').read_text())
    for name, fingerprint in p['source_sha256'].items():
        if digest(ROOT / 'analysis' / name) != fingerprint:
            raise ValueError('Period response implementation changed: ' + name)
    if digest(original.OUT / 'libjointadc.so') != p['unchanged_joint_library_sha256']:
        raise ValueError('The original joint response library changed')
    folder, a, g, adc, domain, train, test = inputs(dataset)
    if digest(folder / 'measured.npz') != p['measurement_sha256'] or digest(folder / 'moments.npz') != p['geometry_sha256']:
        raise ValueError('Period response measurement cache changed')
    spec = p['variants'][variant]
    feature, norm = periodic.features(a, train, dataset, variant, spec['normalization'])
    match = [row for row in p['configs'] if row[0] == family]
    if len(match) != 1:
        raise ValueError('Unknown time/rigidity family')
    _, field, nu = match[0]
    return p, a, g, adc, train, test, feature, field, nu


def verify(dataset):
    destination = OUT / dataset / 'numerical_verification.json'
    if destination.exists():
        raise FileExistsError('Period numeric verification already recorded')
    family = CONFIGS[dataset][0][0]
    p, a, g, adc, train, test, feature, field, nu = context(dataset, family, 'spline')
    seed = p['seeds'][family]
    if digest(seed['path']) != seed['sha256']:
        raise ValueError('Own-dataset numerical seed changed')
    theta = periodic.lift(np.asarray(json.loads(Path(seed['path']).read_text())['model']['theta']), feature.shape[2])
    theta[39:].reshape(4, feature.shape[2])[:, -4:] = np.linspace(-.025, .025, 16).reshape(4, 4)
    selected = np.flatnonzero(train)[np.linspace(0, int(train.sum())-1, 8).astype(int)]
    obj = original.JointADC(a, g, selected, feature, adc, nodes=64, time_field=field, rigidity_nu=nu, threads=1, penalty=.01)
    value, gradient, per = obj.evaluate(theta, events=True)
    numpy_error = np.array([-obj.numpy_event(theta, j)['loglike'] - per[j] for j in range(len(selected))])
    indices = [0, 1, 2, 8, 29, 30, 34, 35] + [39 + layer * feature.shape[2] + column for layer in range(4) for column in range(feature.shape[2]-4, feature.shape[2])]
    differences = []
    for index in indices:
        step = 2e-6; upper = theta.copy(); lower = theta.copy(); upper[index] += step; lower[index] -= step
        finite = (obj.evaluate(upper)[0] - obj.evaluate(lower)[0]) / (2 * step)
        differences.append(abs(finite - gradient[index]))
    if np.max(abs(numpy_error)) > 1e-10 or max(differences) > 3e-5:
        raise ValueError('Period response likelihood or derivative verification failed')
    output = original.reconstruct(theta, g['x'][selected], obj.y, feature[selected])
    injected = np.linspace(-.04, .04, len(selected))
    changed = original.reconstruct(theta, g['x'][selected], obj.y + injected[:, None] * output['direction'], feature[selected])
    response_error = float(np.max(abs(changed['joint_u'] - output['joint_u'] - injected)))
    if response_error > 2e-12:
        raise ValueError('Period response violated combined time/ADC signal equivariance')
    write(destination, dict(complete=True, created_utc=stamp(), dataset=dataset, events=len(selected),
        dimensions=len(theta), gradient_dimensions=len(indices), cpp_numpy_max_nll_error=float(np.max(abs(numpy_error))),
        maximum_central_difference_error=max(differences), signal_injection_max_error=response_error,
        gradient_indices=indices, protocol_sha256=digest(OUT / dataset / 'protocol.json'),
        RICH_or_MC_truth_values_used=False))
    print('PERIOD_ADC_NUMERIC', dataset, len(theta), float(np.max(abs(numpy_error))), max(differences), response_error, flush=True)


def start(dataset, family, variant, number):
    p, a, g, adc, train, test, feature, field, nu = context(dataset, family, variant)
    if number not in range(p['starts']):
        raise ValueError('Invalid period fit initialization')
    numeric = json.loads((OUT / dataset / 'numerical_verification.json').read_text())
    if not numeric['complete'] or numeric['protocol_sha256'] != digest(OUT / dataset / 'protocol.json'):
        raise ValueError('Period response numeric verification required')
    destination = OUT / dataset / (family + '_' + variant); destination.mkdir(exist_ok=True)
    path = destination / ('start_%d.json' % number)
    if path.exists():
        raise FileExistsError('Period fit start already completed')
    seed = p['seeds'][family]
    if digest(seed['path']) != seed['sha256']:
        raise ValueError('Period initialization source changed')
    if number == 0:
        theta = periodic.lift(np.asarray(json.loads(Path(seed['path']).read_text())['model']['theta']), feature.shape[2])
    else:
        theta = original.parameters_initial(feature, adc, train, lambda0=3.)
    obj = original.JointADC(a, g, train, feature, adc, nodes=p['nodes'], time_field=field, rigidity_nu=nu, threads=4, penalty=p['penalty'])
    history = []
    for iteration in range(p['maximum_rounds']):
        result = solve_scaled(obj, theta, p['variants'][variant]['bounds'])
        history.append({key: val for key, val in result.items() if key != 'theta'})
        print('PERIOD_ADC_START', dataset, family, variant, number, iteration, result['penalized_nll'], result['projected_gradient_max'], result['elapsed_seconds'], flush=True)
        if result['success'] and result['projected_gradient_max'] < p['projected_gradient_max']:
            break
        theta = np.asarray(result['theta'])
    result.update(complete=True, created_utc=stamp(), dataset=dataset, family=family, variant=variant, start=number,
        rounds=history, protocol_sha256=digest(OUT / dataset / 'protocol.json'),
        numerically_converged=bool(result['success'] and result['projected_gradient_max'] < p['projected_gradient_max']))
    write(path, result)
    print('PERIOD_ADC_START_SAVED', dataset, family, variant, number, result['numerically_converged'], flush=True)


def freeze(dataset, family, variant):
    p, a, g, adc, train, test, feature, field, nu = context(dataset, family, variant)
    destination = OUT / dataset / (family + '_' + variant)
    if (destination / 'freeze.json').exists():
        raise FileExistsError('Period model already frozen')
    starts = [json.loads((destination / ('start_%d.json' % number)).read_text()) for number in range(p['starts'])]
    if any(not row['complete'] or row['protocol_sha256'] != digest(OUT / dataset / 'protocol.json') for row in starts):
        raise ValueError('Incomplete or changed period fit starts')
    best = min(starts, key=lambda row: row['penalized_nll'])
    if not best['numerically_converged']:
        raise ValueError('Best period fit did not pass convergence; do not silently select a worse successful solution')
    write(destination / 'freeze.json', dict(complete=True, created_utc=stamp(), model=best,
        normalization=p['variants'][variant]['normalization'], starts_sha256={str(number): digest(destination / ('start_%d.json' % number)) for number in range(p['starts'])},
        all_start_convergence=[row['numerically_converged'] for row in starts], protocol_sha256=digest(OUT / dataset / 'protocol.json'),
        RICH_targets_used=False, MC_truth_targets_used=False))
    theta = np.asarray(best['theta'])
    ev = original.JointADC(a, g, test, feature, adc, nodes=p['evaluation_nodes'], time_field=field, rigidity_nu=nu, threads=4, penalty=0.)
    value, _, per = ev.evaluate(theta, events=True)
    ev256 = original.JointADC(a, g, test, feature, adc, nodes=p['nodes'], time_field=field, rigidity_nu=nu, threads=4, penalty=0.)
    value256, _, per256 = ev256.evaluate(theta, events=True)
    result = original.reconstruct(theta, g['x'][test], ev.y, feature[test])
    x = g['x'][test]; ols_weights = x - x.mean(1)[:, None]; ols_weights /= np.sum(ols_weights * x, axis=1)[:, None]
    slopes = dict(ols=(np.sum(ols_weights * a[field][test], axis=1), ols_weights),
        time=(result['time_u'], result['time_weights']), joint=(result['joint_u'], result['weights'][:, :3].dot(timing.D)))
    data = {key: a[key][test] for key in ['run', 'event', 'R']}
    data.update(nll=per, nll256=per256, adc_u=result['adc_u'])
    corrections = {}
    for name, (slope, weights) in slopes.items():
        corrected = transport.solve_center(slope, a['R'][test], weights, g['moment'][test], theta[:2])
        if not corrected['converged'].all():
            raise ValueError('Period mean-loss inverse failed; preserve every event')
        data[name + '_u'] = slope; data[name + '_transport_u'] = corrected['u']; data[name + '_weights'] = weights
        corrections[name] = dict(max_residual=float(np.max(abs(corrected['residual']))), median_correction=float(np.median(corrected['correction'])))
    injected = np.linspace(-.1, .1, len(ev.y))
    changed = original.reconstruct(theta, x, ev.y + injected[:, None] * result['direction'], feature[test])
    error = float(np.max(abs(changed['joint_u'] - result['joint_u'] - injected)))
    if error > 2e-12:
        raise ValueError('Period free-slope coherent injection failed')
    np.savez_compressed(destination / 'predictions_before_reference.npz', **data)
    write(destination / 'measurement_evaluation.json', dict(complete=True, created_utc=stamp(), dataset=dataset, family=family, variant=variant,
        events=len(per), runs=len(np.unique(a['run'][test])), heldout_joint_nll=value,
        quadrature_mean256_minus512=value256-value, quadrature_rms256_minus512=float(np.sqrt(np.mean((per256-per)**2))),
        material_kappa=theta[:2], ADC_lambda=np.exp(theta[35:39]),
        model_information_width_gain_percentiles=np.percentile(1-np.sqrt(result['time_information']/result['information']), [0, 5, 50, 95, 100]),
        time_only_signal_gain_percentiles=np.percentile(result['time_injection_gain'], [0, 5, 50, 95, 100]),
        full_signal_injection_max_error=error, loss_inverse=corrections,
        freeze_sha256=digest(destination / 'freeze.json'), predictions_sha256=digest(destination / 'predictions_before_reference.npz'),
        RICH_targets_used=False, MC_truth_targets_used=False, development_reference_values_scored=False))
    print('PERIOD_ADC_FROZEN', dataset, family, variant, value, 'quadrature', value256-value, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('stage', choices=['prepare', 'verify', 'start', 'freeze'])
    parser.add_argument('--dataset', choices=['mc', 'flight'], required=True)
    parser.add_argument('--family', default='newcalib_gaussianR'); parser.add_argument('--variant', choices=periodic.VARIANTS, default='static')
    parser.add_argument('--start', type=int, default=0); args = parser.parse_args()
    if args.stage in ['prepare', 'verify']:
        globals()[args.stage](args.dataset)
    elif args.stage == 'start':
        start(args.dataset, args.family, args.variant, args.start)
    else:
        freeze(args.dataset, args.family, args.variant)
