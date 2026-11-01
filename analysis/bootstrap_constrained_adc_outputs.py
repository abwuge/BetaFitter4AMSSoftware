#!/usr/bin/env python3
"""Propagate existing measured-run recalibrations through constrained inference."""
import argparse
import json
from pathlib import Path
import numpy as np
import joint_adc_likelihood as joint
import constrained_joint_adc as constrained
import free_timing_transport as transport
from research_constrained_joint_adc import OUT, EXTRA, JOINT, VARIANTS, inputs
from research_flight_conditional import load, write, digest
from research_conditional_mc import ROOT, stamp
from evaluate_additional_mc_20261013 import matrix_metrics, matrix_difference, METRICS, DELTA

BOOT = JOINT / 'bootstrap'
NEW = [variant + suffix for variant in VARIANTS for suffix in ['_u', '_transport_u']]
LABELS = ['native_u', 'time_u', 'joint_u', 'time_transport_u', 'joint_transport_u'] + NEW
PAIRS = [(key, 'native_u') for key in NEW] + [(key, 'time_transport_u' if key.endswith('_transport_u') else 'time_u') for key in NEW]


def plan():
    path = OUT / 'retraining_protocol.json'
    if path.exists():
        raise FileExistsError(path)
    write(path, dict(created_utc=stamp(), training_replicates=128, evaluation_draws=16, seed=911202610,
        source_sha256=digest(Path(__file__)), constrained_sha256=digest(ROOT / 'analysis/constrained_joint_adc.py'),
        original_protocol_sha256=digest(OUT / 'protocol.json'), bootstrap_protocol_sha256=digest(BOOT / 'protocol.json'),
        reference_sha256=digest(EXTRA / 'reference/reference.npz'), branches=NEW, pairs=PAIRS,
        calibration='Apply the two response constraints to each existing genuine run-retrained joint model, using its own recomputed geometry normalization',
        evaluation='Paired whole-run evaluation samples, sixteen draws per training replicate; retain all completed fits',
        scope='Additional C12 development evaluation plus reference-free flight calibration dispersion; no independent C13 or RICH result',
        quantiles=[2.5, 50, 97.5], RICH_read=False))
    print('CONSTRAINED_RETRAINING_PROTOCOL', flush=True)


def outputs(theta, a, g, adc, select, field, norm):
    features, _ = joint.response_features(a, select, norm)
    x = g['x'][select]; yy = joint.observation({field: a[field][select]}, adc[select], field)
    result = {}
    for variant in VARIANTS:
        out = constrained.reconstruct(theta, x, yy, features[select], variant)
        solved = transport.solve_center(out['u'], a['R'][select], out['time_weights'], g['moment'][select], theta[:2])
        if not solved['converged'].all():
            raise ValueError('Constrained bootstrap transport failure')
        result[variant + '_u'] = out['u']; result[variant + '_transport_u'] = solved['u']
    if not all(np.isfinite(v).all() and np.all(v != 0) for v in result.values()):
        raise ValueError('Nonfinite bootstrap prediction')
    return result


def score():
    p = json.loads((OUT / 'retraining_protocol.json').read_text())
    checks = [(Path(__file__), 'source_sha256'), (ROOT / 'analysis/constrained_joint_adc.py', 'constrained_sha256'),
              (OUT / 'protocol.json', 'original_protocol_sha256'), (BOOT / 'protocol.json', 'bootstrap_protocol_sha256'),
              (EXTRA / 'reference/reference.npz', 'reference_sha256')]
    for path, key in checks:
        if digest(path) != p[key]:
            raise ValueError('Constrained bootstrap input changed ' + key)
    if (OUT / 'retraining_evaluation.json').exists():
        raise FileExistsError('Constrained bootstrap already summarized')
    a = load(EXTRA / 'measured.npz'); g = load(EXTRA / 'moments.npz'); select = np.ones(len(a['run']), bool)
    ref = load(EXTRA / 'reference/reference.npz'); original = load(EXTRA / 'predictions_before_truth.npz')
    for key in ['run', 'event']:
        if not np.array_equal(a[key], ref[key]) or not np.array_equal(a[key], original[key]):
            raise ValueError('Additional identities differ')
    valid = ref['valid']; truth = ref['u'][valid]; native = original['native_u'][valid]; runs = ref['run'][valid]
    groups = [np.flatnonzero(runs == r) for r in np.unique(runs)]
    folder, fa, fg, fadc, fdomain, ftrain, ftest = inputs('flight')
    pairs = [(LABELS.index(new), LABELS.index(old)) for new, old in PAIRS]
    samples = {space: dict(train_only=[], joint=[]) for space in ['u', 'beta']}
    fvalues = {key: [] for key in NEW}; evidence = {}; failures = []
    for rep in range(p['training_replicates']):
        path = BOOT / ('%04d.json' % rep); note = json.loads(path.read_text())
        if not note['complete'] or note['replicate'] != rep or note['protocol_sha256'] != p['bootstrap_protocol_sha256'] or digest(path.with_suffix('.npz')) != note['output_sha256']:
            raise ValueError('Invalid joint bootstrap manifest')
        evidence[path.name] = digest(path)
        models = {row['dataset']: row for row in note['models']}
        for dataset, row in models.items():
            if not row['numerically_converged']:
                failures.append(dict(replicate=rep, dataset=dataset))
        own = models['mc']; fp = models['flight']
        new = outputs(np.asarray(own['model']['theta']), a, g, a['adc'], select, 't_native', own['normalization'])
        flight = outputs(np.asarray(fp['model']['theta']), fa, fg, fadc, ftest, 't', fp['normalization'])
        previous = load(path.with_suffix('.npz'))
        for key in ['run', 'event']:
            if not np.array_equal(previous['extra_' + key], a[key]) or not np.array_equal(previous['flight_' + key], fa[key][ftest]):
                raise ValueError('Bootstrap output identity changed')
        values = np.column_stack([native] + [previous['extra_' + key][valid] for key in LABELS[1:5]] + [new[key][valid] for key in NEW])
        rng = np.random.RandomState(p['seed'] + rep)
        draws = [np.concatenate([groups[j] for j in rng.randint(0, len(groups), len(groups))]) for _ in range(p['evaluation_draws'])]
        for space in samples:
            error = values - truth[:, None] if space == 'u' else 1 / values - 1 / truth[:, None]
            met = matrix_metrics(error)
            samples[space]['train_only'].append(np.array([matrix_difference(met[old], met[new]) for new, old in pairs]))
            for ix in draws:
                met = matrix_metrics(error[ix])
                samples[space]['joint'].append(np.array([matrix_difference(met[old], met[new]) for new, old in pairs]))
        for key in NEW:
            fvalues[key].append(flight[key])
        if (rep + 1) % 16 == 0:
            print('CONSTRAINED_RETRAINING', rep + 1, '/', p['training_replicates'], flush=True)
    intervals = {}; arrays = {}
    for space, kinds in samples.items():
        intervals[space] = {}
        for j, (name, old) in enumerate(PAIRS):
            label = name + '_versus_' + old; intervals[space][label] = {}
            for kind, data in kinds.items():
                v = np.asarray(data); arrays[space + '_' + kind] = v
                intervals[space][label][kind] = {k: np.percentile(v[:, j, h], p['quantiles']) for h, k in enumerate(DELTA)}
    flight_summary = {}
    for key, values in fvalues.items():
        sd = np.std(values, axis=0, ddof=1); arrays['flight_' + key + '_std'] = sd
        flight_summary[key] = dict(calibration_std_percentiles=np.percentile(sd, [0, 5, 50, 95, 100]))
    np.savez_compressed(OUT / 'retraining_metrics.npz', **arrays)
    write(OUT / 'retraining_evaluation.json', dict(complete=True, created_utc=stamp(), events=len(truth), runs=len(groups),
        replicates=p['training_replicates'], nested_combinations=p['training_replicates'] * p['evaluation_draws'],
        convergence_failures=failures, numerically_accepted=not failures, intervals=intervals, flight=flight_summary,
        reference_sha256=p['reference_sha256'], protocol_sha256=digest(OUT / 'retraining_protocol.json'),
        metrics_sha256=digest(OUT / 'retraining_metrics.npz'), evidence_sha256=evidence, RICH_read=False))
    for name, old in PAIRS:
        row = intervals['beta'][name + '_versus_' + old]['joint']
        print('CONSTRAINED_RETRAINED_BETA', name, old, 'width', row['w68_improvement'], 'rmse', row['rmse_improvement'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('stage', choices=['plan', 'score']); args = parser.parse_args()
    globals()[args.stage]()
