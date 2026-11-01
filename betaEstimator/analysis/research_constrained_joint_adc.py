#!/usr/bin/env python3
"""Validate ADC noise correction with unit time response on frozen joint models."""
import argparse
import json
from pathlib import Path
import numpy as np
import joint_adc_likelihood as joint
import constrained_joint_adc as constrained
import free_timing_transport as transport
import flight_conditional_likelihood as base
from research_amplitude_noise_20261013 import inputs, CONFIGS
from research_flight_conditional import load, write, digest
from research_conditional_mc import stamp, ROOT
from research_additional_mc_20261013 import OUT as EXTRA
from evaluate_additional_mc_20261013 import matrix_metrics, matrix_difference, METRICS, DELTA

OUT = ROOT / 'results/flight_calibration_research/continuation_20261013/constrained_adc'
JOINT = joint.OUT
VARIANTS = list(constrained.VARIANTS)
FILES = ['constrained_joint_adc.py', 'research_constrained_joint_adc.py', 'joint_adc_likelihood.py', 'free_timing_transport.py']


def check_sources(p):
    for name, value in p['source_sha256'].items():
        if digest(ROOT / 'analysis' / name) != value:
            raise ValueError('Constrained implementation changed ' + name)
    for row in p['models']:
        if digest(row['path']) != row['sha256']:
            raise ValueError('Inherited joint calibration changed')


def plan():
    OUT.mkdir(exist_ok=True)
    if (OUT / 'protocol.json').exists():
        raise FileExistsError('Constrained experiment already declared')
    models = []
    for dataset in ['mc', 'flight']:
        for name, field, nu in CONFIGS[dataset]:
            path = JOINT / dataset / name / 'freeze.json'
            models.append(dict(dataset=dataset, family=name, field=field, path=str(path), sha256=digest(path)))
    write(OUT / 'protocol.json', dict(created_utc=stamp(), models=models, variants=VARIANTS,
        inference='Minimum-variance linear estimator with unit response to Dt speed direction and zero response to fitted ADC speed direction; second variant also nulls common logADC gain',
        objective='Minimize w.T Sigma w with fixed independently calibrated shared noise and mean response',
        proof='For feasible h with C.T h=0, Var(w+h)-Var(w)=h.T Sigma h>=0. This is fixed-model linear-estimator optimality, not measured-detector or W68 global optimality.',
        primary='MC native_studentR5 speed_gain_null_transport_u versus native beta; speed_null and no-transport controls all retained',
        exposure='Developed estimator on already scored C12 runs. Additional-run outcomes for unconstrained estimators are known; these are not new prospective discovery intervals.',
        reference='Use the existing frozen primary-L5 reference only after constrained predictions are saved',
        constants_and_coefficients='No new fit, coefficient tuning, event output cuts or external constants. MC and flight use their own inherited models.',
        branches=['time_u', 'joint_u', 'time_transport_u', 'joint_transport_u'] + [v + suffix for v in VARIANTS for suffix in ['_u', '_transport_u']],
        spaces=['u', 'beta'], metrics=METRICS, deltas=DELTA, replicates=2048, seed=816202610,
        quantiles=[2.5, 50, 97.5], bins=dict(R=[2.5, 3, 3.5, 4.5, 6], beta=[0, .8, .85, .9, .95, 1.01]),
        source_sha256={n: digest(ROOT / 'analysis' / n) for n in FILES},
        test_sha256=digest(ROOT / 'tests/test_constrained_joint_adc.py'),
        extra_input_sha256=digest(EXTRA / 'measured.npz'), extra_geometry_sha256=digest(EXTRA / 'moments.npz'),
        reference_sha256=digest(EXTRA / 'reference/reference.npz'), RICH_read=False))
    print('CONSTRAINED_PROTOCOL', flush=True)


def infer(theta, a, g, adc, features, select, field):
    x = g['x'][select]
    yy = joint.observation({field: a[field][select]}, adc[select], field)
    orig = joint.reconstruct(theta, x, yy, features[select])
    pred = {k: a[k][select] for k in ['run', 'event', 'R']}
    if 'tof_betah' in a:
        pred['native_u'] = 1 / a['tof_betah'][select]
    xc = x - x.mean(1)[:, None]
    pred['ols_u'] = np.sum(xc * a[field][select], axis=1) / np.sum(xc * xc, axis=1)
    rows = {}
    work = [('time', orig['time_u'], orig['time_weights']),
            ('joint', orig['joint_u'], orig['weights'][:, :3].dot(base.D))]
    for variant in VARIANTS:
        result = constrained.reconstruct(theta, x, yy, features[select], variant)
        u = result['u']; w = result['weights']
        work.append((variant, u, result['time_weights']))
        b = np.linspace(-.05, .05, len(u))
        c = np.linspace(-.2, .2, len(u))
        shift = np.column_stack([x.dot(base.D.T) * b[:, None], np.zeros((len(u), 4))])
        time_shift = constrained.reconstruct(theta, x, yy + shift, features[select], variant)['u']
        adc_shift = yy.copy(); adc_shift[:, 3:] += c[:, None] * np.exp(theta[35:39])[None, :]
        adc_u = constrained.reconstruct(theta, x, adc_shift, features[select], variant)['u']
        gain_shift = yy.copy(); gain_shift[:, 3:] += c[:, None]
        gain_u = constrained.reconstruct(theta, x, gain_shift, features[select], variant)['u']
        time_error = float(np.max(abs(time_shift - u - b)))
        adc_error = float(np.max(abs(adc_u - u)))
        gain_response = np.sum(w[:, 3:], axis=1)
        if time_error > 2e-10 or adc_error > 2e-10:
            raise ValueError('Constrained injection failed')
        if variant == 'speed_gain_null' and np.max(abs(gain_u - u)) > 2e-10:
            raise ValueError('Common ADC gain constraint failed')
        rows[variant] = dict(time_injection_max_error=time_error, ADC_speed_injection_max_error=adc_error,
            common_logADC_gain_response_percentiles=np.percentile(gain_response, [0, 5, 50, 95, 100]),
            common_gain_injection_max_change=float(np.max(abs(gain_u - u))),
            constraint_max_error=float(np.max(abs(result['response'] - np.array([1.] + [0.] * (result['response'].shape[1] - 1))))),
            model_width_gain_vs_time_percentiles=np.percentile(result['model_width_gain_vs_time'], [0, 5, 50, 95, 100]),
            gram_condition_percentiles=np.percentile(result['gram_condition'], [0, 50, 100]))
    for name, u, wt in work:
        pred[name + '_u'] = u
        solved = transport.solve_center(u, a['R'][select], wt, g['moment'][select], theta[:2])
        if not solved['converged'].all():
            raise ValueError('Constrained material inverse failed')
        pred[name + '_transport_u'] = solved['u']
        if name in VARIANTS:
            # The material-aware injection follows the same declared forward map.
            b = np.linspace(-.002, .002, len(u))
            m0, m1 = transport.compress(g['moment'][select])
            oldloss = transport.continuous_kernel(solved['u'], a['R'][select], m0, m1)
            newloss = transport.continuous_kernel(solved['u'] + b, a['R'][select], m0, m1)
            delta = np.einsum('nl,nlg,g->n', wt, newloss - oldloss, theta[:2])
            check = transport.solve_center(u + b + delta, a['R'][select], wt, g['moment'][select], theta[:2])
            err = float(np.max(abs(check['u'] - solved['u'] - b)))
            if err > 2e-10:
                raise ValueError('Material-aware speed response failed')
            rows[name]['material_forward_injection_max_error'] = err
            rows[name]['mean_loss_inverse_max_residual'] = float(np.max(abs(solved['residual'])))
    if not all(np.isfinite(v).all() and np.all(v != 0) for k, v in pred.items() if k.endswith('_u')):
        raise ValueError('Invalid prediction; no event-dependent output cuts allowed')
    return pred, rows


def predict():
    p = json.loads((OUT / 'protocol.json').read_text()); check_sources(p)
    dest = OUT / 'output_freeze.json'
    if dest.exists():
        raise FileExistsError('Constrained predictions already frozen')
    if digest(EXTRA / 'measured.npz') != p['extra_input_sha256'] or digest(EXTRA / 'moments.npz') != p['extra_geometry_sha256']:
        raise ValueError('Additional measurements changed')
    evidence = {}; diagnostics = {}
    for dataset in ['mc', 'flight', 'additional_mc']:
        parent = 'mc' if dataset == 'additional_mc' else dataset
        if dataset == 'additional_mc':
            a = load(EXTRA / 'measured.npz'); g = load(EXTRA / 'moments.npz'); adc = a['adc']; test = np.ones(len(a['run']), bool); train = test
        else:
            folder, a, g, adc, domain, train, test = inputs(dataset)
        cp = json.loads((JOINT / parent / 'protocol.json').read_text())
        features, _ = joint.response_features(a, train, cp['feature_normalization'])
        for row in [r for r in p['models'] if r['dataset'] == parent]:
            theta = np.array(json.loads(Path(row['path']).read_text())['model']['theta'])
            pred, diag = infer(theta, a, g, adc, features, test, row['field'])
            label = dataset + '_' + row['family']; path = OUT / (label + '.npz')
            np.savez_compressed(path, **pred)
            evidence[path.name] = digest(path); diagnostics[label] = diag
            print('CONSTRAINED_FROZEN', label, len(pred['run']), diag, flush=True)
    write(dest, dict(complete=True, created_utc=stamp(), protocol_sha256=digest(OUT / 'protocol.json'),
        predictions_sha256=evidence, injection_diagnostics=diagnostics, MC_truth_read=False, RICH_read=False))


def score():
    p = json.loads((OUT / 'protocol.json').read_text()); check_sources(p)
    gate = json.loads((OUT / 'output_freeze.json').read_text())
    if gate['protocol_sha256'] != digest(OUT / 'protocol.json'):
        raise ValueError('Constrained output protocol changed')
    for name, value in gate['predictions_sha256'].items():
        if digest(OUT / name) != value:
            raise ValueError('Constrained predictions changed')
    dest = OUT / 'evaluation.json'
    if dest.exists():
        raise FileExistsError('Constrained score already exists')
    if digest(EXTRA / 'reference/reference.npz') != p['reference_sha256']:
        raise ValueError('Primary reference changed')
    ref = load(EXTRA / 'reference/reference.npz'); valid = ref['valid']; truth = ref['u'][valid]; runs = ref['run'][valid]
    arrays = {}; labels = []; first = None
    for family in ['native_gaussianR', 'native_studentR5']:
        pred = load(OUT / ('additional_mc_' + family + '.npz'))
        for key in ['run', 'event']:
            if not np.array_equal(pred[key], ref[key]):
                raise ValueError('Constrained scoring identity mismatch')
        if first is None:
            first = pred; labels.append('native_u'); arrays['native_u'] = pred['native_u'][valid]
        for branch in p['branches']:
            name = family + '_' + branch; labels.append(name); arrays[name] = pred[branch][valid]
    values = np.column_stack([arrays[k] for k in labels]); groups = [np.flatnonzero(runs == r) for r in np.unique(runs)]
    pairs = [(k, 'native_u') for k in labels[1:]]
    for family in ['native_gaussianR', 'native_studentR5']:
        for variant in VARIANTS:
            for suffix in ['_u', '_transport_u']:
                pairs.append((family + '_' + variant + suffix, family + '_time' + suffix))
                pairs.append((family + '_' + variant + suffix, family + '_joint' + suffix))
    point = {}; samples = {}; rng = np.random.RandomState(p['seed'])
    errors = {space: values - truth[:, None] if space == 'u' else 1 / values - 1 / truth[:, None] for space in p['spaces']}
    for space in errors:
        point[space] = matrix_metrics(errors[space]); samples[space] = np.zeros((p['replicates'], len(labels), len(METRICS)))
    for rep in range(p['replicates']):
        index = np.concatenate([groups[j] for j in rng.randint(0, len(groups), len(groups))])
        for space in errors:
            samples[space][rep] = matrix_metrics(errors[space][index])
        if (rep + 1) % 512 == 0:
            print('CONSTRAINED_EVALUATION_BOOTSTRAP', rep + 1, flush=True)
    table = {}; bins = []
    for space in p['spaces']:
        comparisons = {}
        for name, baseline in pairs:
            j = labels.index(name); b = labels.index(baseline)
            delta = matrix_difference(point[space][b], point[space][j]); rep = matrix_difference(samples[space][:, b], samples[space][:, j])
            comparisons[name + '_versus_' + baseline] = dict(point=dict(zip(DELTA, delta)),
                interval={k: np.percentile(rep[:, h], p['quantiles']) for h, k in enumerate(DELTA)})
        table[space] = dict(metrics={k: dict(zip(METRICS, point[space][j])) for j, k in enumerate(labels)}, comparisons=comparisons)
        for axis, edges in p['bins'].items():
            coordinate = first['R'][valid] if axis == 'R' else 1 / truth
            for low, high in zip(edges[:-1], edges[1:]):
                take = (coordinate >= low) & (coordinate < high)
                if take.sum() < 20:
                    continue
                met = matrix_metrics(errors[space][take]); dif = matrix_difference(met[0][None, :], met)
                bins.append(dict(space=space, axis=axis, limits=[low, high], events=int(take.sum()),
                    metrics={k: dict(zip(METRICS, met[j])) for j, k in enumerate(labels)},
                    versus_native={k: dict(zip(DELTA, dif[j])) for j, k in enumerate(labels)}))
    np.savez_compressed(OUT / 'evaluation_bootstrap.npz', **samples)
    write(dest, dict(complete=True, created_utc=stamp(), events=len(truth), runs=len(groups), table=table, bins=bins,
        reference_sha256=p['reference_sha256'], output_freeze_sha256=digest(OUT / 'output_freeze.json'),
        bootstrap_sha256=digest(OUT / 'evaluation_bootstrap.npz'), labels=labels, RICH_read=False))
    for space in p['spaces']:
        for variant in VARIANTS:
            for suffix in ['_u', '_transport_u']:
                name = 'native_studentR5_' + variant + suffix
                row = table[space]['comparisons'][name + '_versus_native_u']
                print('CONSTRAINED_RESULT', space, name, 'point', row['point'], 'widthCI', row['interval']['w68_improvement'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('stage', choices=['plan', 'predict', 'score'])
    args = parser.parse_args(); globals()[args.stage]()
