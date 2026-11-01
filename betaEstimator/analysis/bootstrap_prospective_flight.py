#!/usr/bin/env python3
"""Replay the saved flight run refits on the fixed prospective RICH cohort."""
import argparse
import json
from pathlib import Path
import numpy as np
from prepare_prospective_flight import ROOT, OUT
from reference_prospective_flight import REF, gate
from predict_prospective_flight import infer_family, JOINT
from research_flight_conditional import load, write, digest
from research_conditional_mc import stamp
from evaluate_additional_mc_20261013 import matrix_metrics, matrix_difference, METRICS, DELTA
import free_timing_transport as center
import continuous_rich_transport as forward

SCORE = OUT / 'scoring'
DEST = SCORE / 'recalibration'
BOOT = JOINT.parent / 'bootstrap'
QUAD = JOINT.parent / 'quadrature_refit/flight/newcalib_gaussianR'
BRANCHES = ['ols', 'time', 'joint', 'speed_null', 'speed_gain_null']
PAIRS = [(branch, 'ols') for branch in BRANCHES if branch != 'ols'] + [(branch, 'time') for branch in ['joint', 'speed_null', 'speed_gain_null']]
SOURCES = ['bootstrap_prospective_flight.py', 'predict_prospective_flight.py', 'joint_adc_likelihood.py',
           'constrained_joint_adc.py', 'free_timing_transport.py', 'continuous_rich_transport.py', 'evaluate_additional_mc_20261013.py']


def plan():
    DEST.mkdir(parents=True, exist_ok=True)
    if (DEST / 'protocol.json').exists():
        raise FileExistsError('Prospective recalibration protocol already exists')
    if (REF / 'read_started.json').exists():
        raise ValueError('Declare recalibration evaluation before opening the prospective RICH outcomes')
    records = {}
    failures = []
    protocol_hash = digest(BOOT / 'protocol.json')
    for rep in range(128):
        path = BOOT / ('%04d.json' % rep)
        note = json.loads(path.read_text())
        rows = [r for r in note['models'] if r['dataset'] == 'flight']
        if not note['complete'] or note['replicate'] != rep or note['protocol_sha256'] != protocol_hash or len(rows) != 1:
            raise ValueError('Saved flight recalibration is incomplete')
        if rows[0]['family'] != 'newcalib_gaussianR':
            raise ValueError('Wrong flight recalibration family')
        if not rows[0]['numerically_converged']:
            failures.append(rep)
        records[path.name] = digest(path)
    point = {'128': dict(path=str(JOINT / 'newcalib_gaussianR/freeze.json'), sha256=digest(JOINT / 'newcalib_gaussianR/freeze.json'))}
    for nodes in [256, 512]:
        path = QUAD / str(nodes) / 'fit.json'
        point[str(nodes)] = dict(path=str(path), sha256=digest(path))
    p = json.loads((OUT / 'prediction_protocol.json').read_text())
    write(DEST / 'protocol.json', dict(created_utc=stamp(), training_replicates=128, evaluation_draws=16,
        seed=1003202611, quantiles=[2.5, 50, 97.5], two_primary_adjusted_quantiles=[1.25, 50, 98.75],
        family='newcalib_gaussianR', branches=BRANCHES, pairs=PAIRS, primary_branches=['joint', 'speed_gain_null'],
        input_sources_sha256={name: digest(ROOT / 'analysis' / name) for name in SOURCES},
        saved_training_protocol_sha256=protocol_hash, saved_refit_sha256=records, inherited_numerical_failures=failures,
        prediction_protocol_sha256=digest(OUT / 'prediction_protocol.json'), original_normalization_sha256=digest(JOINT / 'protocol.json'),
        point_models=point, primary_endpoint=p['primary_endpoint'],
        reconstruction='Replay each saved flight-only run refit with its own geometry normalization; fitted two-group L5 correction and unit forward NaF material for every paired candidate and baseline',
        point_comparison='Report original 128-node, unresampled 256-node, and unresampled 512-node fits separately; never select among them by the new RICH score',
        intervals='Training-only spread and paired training-plus-evaluation-run spread; sixteen independent evaluation-run draws per saved calibration',
        coverage='Preserve every reference event and every saved recalibration. If a branch fails full cohort coverage in any replicate, do not report a filtered all-replicate confidence interval for that branch.',
        uncertainty_limit='Does not resample upstream TOF calibration, material-map systematics, selection-model errors, common RICH-reference errors, or the entire historical model search',
        no_new_training=True, new_RICH_values_read=False))
    (DEST / 'parts').mkdir(exist_ok=True)
    print('PROSPECTIVE_RECALIBRATION_PROTOCOL', 128, 'saved flight refits', 'numerical_failures', failures, flush=True)


def check():
    gate()
    p = json.loads((DEST / 'protocol.json').read_text())
    for name, fingerprint in p['input_sources_sha256'].items():
        if digest(ROOT / 'analysis' / name) != fingerprint:
            raise ValueError('Prospective recalibration source changed: ' + name)
    if digest(JOINT / 'protocol.json') != p['original_normalization_sha256']:
        raise ValueError('Flight normalization changed')
    if digest(BOOT / 'protocol.json') != p['saved_training_protocol_sha256']:
        raise ValueError('Saved calibration protocol changed')
    audit = json.loads((REF / 'audit.json').read_text())
    if not audit['complete'] or digest(REF / 'reference.npz') != audit['reference_sha256']:
        raise ValueError('New native RICH reference is incomplete')
    return p, audit


def inputs():
    a = load(OUT / 'measured.npz')
    g = load(OUT / 'moments.npz')
    rg = load(OUT / 'rich_geometry.npz')
    ref = load(REF / 'reference.npz')
    for key in ['run', 'event']:
        if not np.array_equal(a[key], g[key]) or not np.array_equal(a[key], rg[key]) or not np.array_equal(a[key], ref[key]):
            raise ValueError('Prospective recalibration identity mismatch')
    mask = ref['valid'].astype(bool)
    if mask.sum() < 20:
        raise ValueError('Insufficient prospective reference events')
    return ({k: v[mask] for k, v in a.items()}, {k: v[mask] for k, v in g.items()},
            forward.subset_geometry(rg, mask), {k: v[mask] for k, v in ref.items()}, mask)


def outputs(theta, norm, a, g, rg):
    branches, responses = infer_family(theta, a, g, norm, 't')
    values = []
    coverage = []
    for branch in BRANCHES:
        slope, weights = branches[branch]
        result = center.solve_center(slope, a['R'], weights, g['moment'], theta[:2])
        accepted = result['converged'].copy()
        if not accepted.all():
            values.append(np.full(len(slope), np.nan))
            coverage.append(dict(branch=branch, complete=False, inverse_failures=int((~accepted).sum())))
            continue
        propagated = forward.propagate(result['u'], a['R'], rg)
        u = propagated['u']
        accepted &= propagated['valid'] & np.isfinite(u) & (u != 0)
        values.append(u)
        coverage.append(dict(branch=branch, complete=bool(accepted.all()), invalid_events=int((~accepted).sum())))
    return np.column_stack(values), coverage, responses


def point():
    p, audit = check()
    if (DEST / 'point.json').exists():
        raise FileExistsError('Prospective quadrature point comparison already exists')
    a, g, rg, ref, mask = inputs()
    original = load(OUT / 'predictions_before_reference.npz')
    norm = json.loads((JOINT / 'protocol.json').read_text())['feature_normalization']
    records = {}
    arrays = dict(run=ref['run'], event=ref['event'])
    original_match = {}
    for nodes in ['128', '256', '512']:
        row = p['point_models'][nodes]
        if digest(row['path']) != row['sha256']:
            raise ValueError('Frozen quadrature point model changed')
        fit = json.loads(Path(row['path']).read_text())
        fit = fit['model'] if nodes == '128' else fit
        values, coverage, response = outputs(np.asarray(fit['theta']), norm, a, g, rg)
        arrays['nodes_' + nodes + '_u'] = values
        complete = [c['complete'] for c in coverage]
        valid_branches = [j for j, good in enumerate(complete) if good]
        tables = {}
        for space in ['u', 'beta']:
            value = values[:, valid_branches]
            error = value - 1 / ref['beta'][:, None] if space == 'u' else 1 / value - ref['beta'][:, None]
            met = matrix_metrics(error)
            tables[space] = {BRANCHES[j]: dict(zip(METRICS, met[k])) for k, j in enumerate(valid_branches)}
        if nodes == '128':
            for j, branch in enumerate(BRANCHES):
                expected = original['newcalib_gaussianR_' + branch + '_fit_RICH_s1_u'][mask]
                diff = float(np.nanmax(abs(expected - values[:, j])))
                original_match[branch] = diff
                if diff > 2e-11:
                    raise ValueError('Original prospective reconstruction did not replay exactly')
        records[nodes] = dict(coverage=coverage, metrics=tables, response=response,
            original_optimizer_success=fit.get('success'), original_projected_gradient=fit.get('projected_gradient_max'),
            source_sha256=row['sha256'])
    for nodes in ['256', '512']:
        diff = arrays['nodes_' + nodes + '_u'] - arrays['nodes_128_u']
        records[nodes]['u_difference_from_128'] = {branch: dict(rms=float(np.sqrt(np.nanmean(diff[:, j] ** 2))),
            max_abs=float(np.nanmax(abs(diff[:, j])))) for j, branch in enumerate(BRANCHES)}
    np.savez_compressed(DEST / 'point_predictions.npz', **arrays)
    write(DEST / 'point.json', dict(complete=True, created_utc=stamp(), events=len(ref['run']), records=records,
        original_prediction_replay_max_abs=original_match, reference_sha256=audit['reference_sha256'],
        protocol_sha256=digest(DEST / 'protocol.json'), output_sha256=digest(DEST / 'point_predictions.npz'),
        selection_by_new_score=False))
    print('PROSPECTIVE_QUADRATURE_POINTS', json.dumps({k: records[k]['metrics']['beta'] for k in records}), flush=True)


def replicate(number):
    p, audit = check()
    if number < 0 or number >= p['training_replicates']:
        raise ValueError('Invalid saved calibration index')
    path = DEST / 'parts' / ('%04d.json' % number)
    if path.exists():
        old = json.loads(path.read_text())
        if old['complete'] and old['protocol_sha256'] == digest(DEST / 'protocol.json') and digest(path.with_suffix('.npz')) == old['output_sha256']:
            print('PROSPECTIVE_RECALIBRATION_ALREADY_COMPLETE', number, flush=True)
            return
        raise ValueError('Existing recalibration result failed integrity')
    if path.with_suffix('.npz').exists():
        raise FileExistsError('Unrecorded prospective recalibration output; inspect it before retrying')
    saved = BOOT / ('%04d.json' % number)
    if digest(saved) != p['saved_refit_sha256'][saved.name]:
        raise ValueError('Saved flight calibration changed')
    note = json.loads(saved.read_text())
    own = next(r for r in note['models'] if r['dataset'] == 'flight')
    a, g, rg, ref, mask = inputs()
    values, coverage, response = outputs(np.asarray(own['model']['theta']), own['normalization'], a, g, rg)
    accepted = [j for j, row in enumerate(coverage) if row['complete']]
    arrays = dict(run=ref['run'], event=ref['event'], u=values)
    groups = [np.flatnonzero(ref['run'] == run) for run in np.unique(ref['run'])]
    if len(groups) < 2:
        raise ValueError('At least two reference runs required')
    rng = np.random.RandomState(p['seed'] + number)
    draws = [np.concatenate([groups[j] for j in rng.randint(0, len(groups), len(groups))]) for _ in range(p['evaluation_draws'])]
    for space in ['u', 'beta']:
        error = values[:, accepted] - 1 / ref['beta'][:, None] if space == 'u' else 1 / values[:, accepted] - ref['beta'][:, None]
        met = np.full((len(BRANCHES), len(METRICS)), np.nan)
        samples = np.full((p['evaluation_draws'], len(BRANCHES), len(METRICS)), np.nan)
        met[accepted] = matrix_metrics(error)
        for j, ix in enumerate(draws):
            samples[j, accepted] = matrix_metrics(error[ix])
        arrays[space + '_train_metrics'] = met
        arrays[space + '_joint_metrics'] = samples
    np.savez_compressed(path.with_suffix('.npz'), **arrays)
    write(path, dict(complete=True, created_utc=stamp(), replicate=number, events=len(ref['run']), runs=len(groups),
        coverage=coverage, response=response, numerical_training_converged=own['numerically_converged'],
        saved_model_sha256=p['saved_refit_sha256'][saved.name], reference_sha256=audit['reference_sha256'],
        protocol_sha256=digest(DEST / 'protocol.json'), output_sha256=digest(path.with_suffix('.npz'))))
    print('PROSPECTIVE_RECALIBRATION_COMPLETE', number, flush=True)


def collect():
    p, audit = check()
    if (DEST / 'evaluation.json').exists():
        raise FileExistsError('Prospective recalibration intervals already saved')
    samples = {space: dict(train_only=[], joint=[]) for space in ['u', 'beta']}
    evidence = {}
    bad_coverage = []
    failed_fits = []
    fingerprint = digest(DEST / 'protocol.json')
    for rep in range(p['training_replicates']):
        path = DEST / 'parts' / ('%04d.json' % rep)
        note = json.loads(path.read_text())
        if not note['complete'] or note['protocol_sha256'] != fingerprint or note['reference_sha256'] != audit['reference_sha256'] or digest(path.with_suffix('.npz')) != note['output_sha256']:
            raise ValueError('Incomplete or changed prospective calibration replay')
        evidence[path.name] = digest(path)
        if not note['numerical_training_converged']:
            failed_fits.append(rep)
        for row in note['coverage']:
            if not row['complete']:
                bad_coverage.append(dict(replicate=rep, **row))
        data = load(path.with_suffix('.npz'))
        for space in samples:
            samples[space]['train_only'].append(data[space + '_train_metrics'])
            samples[space]['joint'].extend(data[space + '_joint_metrics'])
    intervals = {}
    arrays = {}
    for space, kinds in samples.items():
        intervals[space] = {}
        for kind, values in kinds.items():
            values = np.asarray(values)
            arrays[space + '_' + kind] = values
            rows = {}
            for new, old in PAIRS:
                newer, older = values[:, BRANCHES.index(new)], values[:, BRANCHES.index(old)]
                complete = bool(np.isfinite(newer).all() and np.isfinite(older).all())
                label = new + '_versus_' + old
                if not complete:
                    rows[label] = dict(complete=False, reason='One or more fixed-cohort calibration predictions invalid; no replicate or event was silently dropped')
                    continue
                change = matrix_difference(older, newer)
                rows[label] = dict(complete=True,
                    interval={key: np.percentile(change[:, j], p['quantiles']) for j, key in enumerate(DELTA)},
                    two_primary_adjusted_interval={key: np.percentile(change[:, j], p['two_primary_adjusted_quantiles']) for j, key in enumerate(DELTA)})
            intervals[space][kind] = rows
    np.savez_compressed(DEST / 'bootstrap_metrics.npz', **arrays)
    write(DEST / 'evaluation.json', dict(complete=True, created_utc=stamp(), replicates=p['training_replicates'],
        nested_combinations=p['training_replicates'] * p['evaluation_draws'], numerical_training_failures=failed_fits,
        coverage_failures=bad_coverage, intervals=intervals, parts_sha256=evidence,
        reference_sha256=audit['reference_sha256'], protocol_sha256=fingerprint,
        metrics_sha256=digest(DEST / 'bootstrap_metrics.npz'), new_model_training=False))
    for branch in p['primary_branches']:
        print('PROSPECTIVE_RETRAINED_PRIMARY', branch, json.dumps(intervals['beta']['joint'][branch + '_versus_ols'],
            default=lambda v: v.tolist() if hasattr(v, 'tolist') else v), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['plan', 'point', 'replicate', 'collect'])
    parser.add_argument('--replicate', type=int)
    args = parser.parse_args()
    replicate(args.replicate) if args.stage == 'replicate' else globals()[args.stage]()
