#!/usr/bin/env python3
"""Summarize preserved run-retrained ADC/time estimators and numerical sensitivity."""
import argparse
import json
from pathlib import Path
import numpy as np
from bootstrap_joint_adc_20261013 import OUT as BOOT, JOINT, EXTRA, prediction, inputs
from research_flight_conditional import load, write, digest
from research_conditional_mc import stamp
from evaluate_additional_mc_20261013 import matrix_metrics, matrix_difference, METRICS, DELTA
import joint_adc_likelihood as joint

SCORE = JOINT / 'scoring'
BRANCHES = ['time_u', 'adc_u', 'joint_u', 'time_transport_u', 'joint_transport_u']
PAIRS = [(k, 'native_u') for k in BRANCHES] + [('joint_u', 'time_u'), ('joint_transport_u', 'time_transport_u')]


def prepare():
    dest = SCORE / 'retraining_protocol.json'
    if dest.exists():
        raise FileExistsError(dest)
    p = json.loads((BOOT / 'protocol.json').read_text())
    write(dest, dict(created_utc=stamp(), replicates=p['replicates'], evaluation_draws_per_training=16,
        seed=613202610, spaces=['u', 'beta'], metrics=METRICS, deltas=DELTA, branches=BRANCHES, pairs=PAIRS,
        primary='MC native_studentR5_joint_transport_u versus native_u, inherited from the pre-existing scoring protocol',
        exposure='Point results on additional C12 runs were already inspected; this protocol fixes only the missing uncertainty summary, not a new blind discovery test',
        training='128 saved whole-training-run recalibrations at 256 quadrature nodes, including recalculated normalization and two optimization starts',
        evaluation='16 independent paired whole-evaluation-run draws per recalibration; 2048 nested combinations, 128 distinct training draws',
        quantiles=[2.5, 50, 97.5], failures='Keep every completed replicate; any convergence failure prevents numerical acceptance',
        flight='Prediction dispersion on the held-out measured flight events; no RICH or beta reference',
        sensitivity='Compute full-training 256-node output and compare it with the original 128-node frozen output; do not relabel original estimates',
        bootstrap_protocol_sha256=digest(BOOT / 'protocol.json'),
        reference_sha256=digest(EXTRA / 'reference/reference.npz'),
        extra_prediction_sha256=digest(EXTRA / 'predictions_before_truth.npz'),
        source_sha256=digest(Path(__file__)), RICH_read=False))
    print('RETRAINING_SUMMARY_PROTOCOL', flush=True)


def score():
    protocol_path = SCORE / 'retraining_protocol.json'
    p = json.loads(protocol_path.read_text())
    if digest(Path(__file__)) != p['source_sha256']:
        raise ValueError('Summary source changed')
    if digest(BOOT / 'protocol.json') != p['bootstrap_protocol_sha256']:
        raise ValueError('Bootstrap protocol changed')
    if digest(EXTRA / 'reference/reference.npz') != p['reference_sha256']:
        raise ValueError('Reference changed')
    if digest(EXTRA / 'predictions_before_truth.npz') != p['extra_prediction_sha256']:
        raise ValueError('Original predictions changed')
    dest = SCORE / 'retraining_evaluation.json'
    if dest.exists():
        raise FileExistsError(dest)
    bp = json.loads((BOOT / 'protocol.json').read_text())
    notes = []
    failures = []
    evidence = {}
    for rep in range(p['replicates']):
        path = BOOT / ('%04d.json' % rep)
        note = json.loads(path.read_text())
        if not note['complete'] or note['replicate'] != rep or note['protocol_sha256'] != p['bootstrap_protocol_sha256']:
            raise ValueError('Bad bootstrap manifest %d' % rep)
        if digest(path.with_suffix('.npz')) != note['output_sha256']:
            raise ValueError('Bootstrap prediction hash mismatch %d' % rep)
        if sorted(m['dataset'] for m in note['models']) != ['flight', 'mc']:
            raise ValueError('Missing bootstrap model')
        for row in note['models']:
            model = row['model']
            converged = model['success'] and model['projected_gradient_max'] < 1e-4
            if bool(converged) != bool(row['numerically_converged']):
                raise ValueError('Inconsistent convergence status')
            if not converged:
                failures.append(dict(replicate=rep, dataset=row['dataset'], pg=model['projected_gradient_max']))
        evidence[path.name] = digest(path)
        notes.append(note)
    original = load(EXTRA / 'predictions_before_truth.npz')
    ref = load(EXTRA / 'reference/reference.npz')
    data = load(EXTRA / 'measured.npz')
    geometry = load(EXTRA / 'moments.npz')
    for key in ['run', 'event']:
        if not np.array_equal(original[key], ref[key]) or not np.array_equal(original[key], data[key]):
            raise ValueError('Additional evaluation identity mismatch')
    valid = ref['valid']
    truth = ref['u'][valid]
    native = original['native_u'][valid]
    runs = ref['run'][valid]
    groups = [np.flatnonzero(runs == r) for r in np.unique(runs)]
    labels = ['native_u'] + BRANCHES
    pairs = [(labels.index(a), labels.index(b)) for a, b in PAIRS]
    samples = {space: dict(train_only=[], joint=[]) for space in p['spaces']}
    flight_outputs = {k: [] for k in BRANCHES}
    parameter_values = {'mc': [], 'flight': []}
    for rep, note in enumerate(notes):
        pred = load(BOOT / ('%04d.npz' % rep))
        for key in ['run', 'event']:
            if not np.array_equal(pred['extra_' + key], original[key]):
                raise ValueError('Bootstrap additional identities changed')
        values = np.column_stack([native] + [pred['extra_' + key][valid] for key in BRANCHES])
        if not np.isfinite(values).all() or np.any(values == 0):
            raise ValueError('Invalid bootstrap output; event removal is forbidden')
        rng = np.random.RandomState(p['seed'] + rep)
        draws = [np.concatenate([groups[j] for j in rng.randint(0, len(groups), len(groups))])
                 for _ in range(p['evaluation_draws_per_training'])]
        for space in p['spaces']:
            errors = values - truth[:, None] if space == 'u' else 1 / values - 1 / truth[:, None]
            point = matrix_metrics(errors)
            samples[space]['train_only'].append(np.array([matrix_difference(point[b], point[a]) for a, b in pairs]))
            for index in draws:
                metrics = matrix_metrics(errors[index])
                samples[space]['joint'].append(np.array([matrix_difference(metrics[b], metrics[a]) for a, b in pairs]))
        for key in BRANCHES:
            flight_outputs[key].append(pred['flight_' + key])
        for row in note['models']:
            parameter_values[row['dataset']].append(row['model']['theta'])
        if (rep + 1) % 32 == 0:
            print('RETRAINING_SUMMARY', rep + 1, '/', p['replicates'], flush=True)
    intervals = {}
    arrays = {}
    for space in p['spaces']:
        intervals[space] = {}
        for j, (a, b) in enumerate(PAIRS):
            label = a + '_versus_' + b
            intervals[space][label] = {}
            for kind, raw in samples[space].items():
                raw = np.asarray(raw)
                arrays[space + '__' + kind] = raw
                intervals[space][label][kind] = {name: np.percentile(raw[:, j, k], p['quantiles']) for k, name in enumerate(DELTA)}
    sensitivity = {}
    for dataset, family in [('mc', 'native_studentR5'), ('flight', 'newcalib_gaussianR')]:
        directory = JOINT / dataset
        cp = json.loads((directory / 'protocol.json').read_text())
        theta = np.array(json.loads((JOINT / 'quadrature_refit' / dataset / family / '256/fit.json').read_text())['theta'])
        if dataset == 'mc':
            selected = np.ones(len(data['run']), bool)
            features, _ = joint.response_features(data, selected, cp['feature_normalization'])
            dense = prediction(theta, data, geometry, data['adc'], features, selected, 't_native')
            old = {key: original[family + '_' + key] for key in BRANCHES}
        else:
            folder, own, geo, adc, domain, train, test = inputs('flight')
            features, _ = joint.response_features(own, train, cp['feature_normalization'])
            dense = prediction(theta, own, geo, adc, features, test, 't')
            old = load(directory / family / 'predictions_before_truth.npz')
            for key in ['run', 'event']:
                if not np.array_equal(dense[key], pred['flight_' + key]):
                    raise ValueError('Flight bootstrap identities changed')
        sensitivity[dataset] = {}
        for key in BRANCHES:
            change = dense[key] - old[key]
            row = dict(rms=float(np.sqrt(np.mean(change * change))), max_abs=float(np.max(abs(change))), median=float(np.median(change)))
            if dataset == 'mc':
                v = np.column_stack([native, dense[key][valid]])
                row['dense_point'] = {}
                for space in p['spaces']:
                    err = v - truth[:, None] if space == 'u' else 1 / v - 1 / truth[:, None]
                    met = matrix_metrics(err)
                    row['dense_point'][space] = dict(metrics=dict(zip(METRICS, met[1])), versus_native=dict(zip(DELTA, matrix_difference(met[0], met[1]))))
            else:
                ensemble = np.asarray(flight_outputs[key])
                sd = np.std(ensemble, axis=0, ddof=1)
                drift = ensemble.mean(0) - dense[key]
                row['recalibration_event_std_percentiles'] = np.percentile(sd, [0, 5, 50, 95, 100])
                row['ensemble_mean_minus_dense_rms'] = float(np.sqrt(np.mean(drift * drift)))
                arrays['flight_' + key + '_calibration_std'] = sd
            sensitivity[dataset][key] = row
    parameters = {}
    for dataset, rows in parameter_values.items():
        theta = np.asarray(rows)
        parameters[dataset] = dict(material_percentiles=np.percentile(theta[:, :2], p['quantiles'], axis=0),
            ADC_slope_percentiles=np.percentile(np.exp(theta[:, 35:39]), p['quantiles'], axis=0))
    np.savez_compressed(SCORE / 'retraining_metrics.npz', **arrays)
    write(dest, dict(complete=True, created_utc=stamp(), additional_events=len(truth), additional_runs=len(groups),
        replicates=len(notes), nested_combinations=len(notes) * p['evaluation_draws_per_training'],
        convergence_failures=failures, numerically_accepted=not failures, intervals=intervals,
        quadrature_sensitivity=sensitivity, parameters=parameters, evidence_sha256=evidence,
        protocol_sha256=digest(protocol_path), metrics_sha256=digest(SCORE / 'retraining_metrics.npz'),
        reference_sha256=p['reference_sha256'], RICH_read=False))
    for space in p['spaces']:
        for a, b in PAIRS:
            row = intervals[space][a + '_versus_' + b]['joint']
            print('RETRAINING_INTERVAL', space, a, b, 'W68', row['w68_improvement'], 'RMSE', row['rmse_improvement'], flush=True)
    print('FLIGHT_CALIBRATION_DISPERSION', sensitivity['flight']['joint_transport_u'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['prepare', 'score'])
    args = parser.parse_args()
    prepare() if args.stage == 'prepare' else score()
