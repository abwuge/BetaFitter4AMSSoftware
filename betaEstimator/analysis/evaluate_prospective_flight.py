#!/usr/bin/env python3
"""Score the prospectively frozen flight estimators on one fixed RICH cohort."""
import argparse
import json
from datetime import datetime
from pathlib import Path
import numpy as np
from prepare_prospective_flight import ROOT, OUT
from reference_prospective_flight import REF, gate
from research_flight_conditional import load, write, digest
from research_conditional_mc import stamp
from evaluate_additional_mc_20261013 import matrix_metrics, matrix_difference, METRICS, DELTA

SCORE = OUT / 'scoring'
FAMILIES = ['native_gaussianR', 'newcalib_gaussianR', 'native_studentR5']
BRANCHES = ['ols', 'time', 'joint', 'speed_null', 'speed_gain_null']
NOMINAL = [family + '_' + branch + '_fit_RICH_s1_u' for family in FAMILIES for branch in BRANCHES]
CONTROLS = ['native_BetaH_u', 'constant_one_u']


def plan():
    SCORE.mkdir(exist_ok=True)
    if (SCORE / 'protocol.json').exists():
        raise FileExistsError('Prospective score protocol already declared')
    if (REF / 'read_started.json').exists():
        raise ValueError('Declare the score implementation before opening the new RICH outcomes')
    p = json.loads((OUT / 'prediction_protocol.json').read_text())
    pairs = []
    for family in FAMILIES:
        baseline = family + '_ols_fit_RICH_s1_u'
        time = family + '_time_fit_RICH_s1_u'
        pairs.extend([[family + '_' + branch + '_fit_RICH_s1_u', baseline] for branch in BRANCHES if branch != 'ols'])
        pairs.extend([[family + '_' + branch + '_fit_RICH_s1_u', time] for branch in ['joint', 'speed_null', 'speed_gain_null']])
    write(SCORE / 'protocol.json', dict(created_utc=stamp(), source_sha256=digest(Path(__file__)),
        metric_source_sha256=digest(ROOT / 'analysis/evaluate_additional_mc_20261013.py'),
        prediction_protocol_sha256=digest(OUT / 'prediction_protocol.json'),
        reference_implementation_sha256=digest(ROOT / 'analysis/reference_prospective_flight.py'),
        primary_comparisons=p['primary_comparisons'], primary_endpoint=p['primary_endpoint'],
        nominal_outputs=NOMINAL + CONTROLS, pairs=pairs, metrics=METRICS, differences=DELTA,
        point_evaluation='All frozen RICH material variants and the two declared untransported controls; identical native-qualified RICH cohort for every valid branch',
        bootstrap=p['evaluation_bootstrap'],
        bootstrap_scope='Nominal fitted-material, unit-forward-scale outputs and their predeclared paired comparisons; all material sensitivities are reported as point diagnostics without choosing the best variant',
        bins=p['bins'], minimum_bin_events=20,
        coverage='A candidate with any nonfinite, zero or invalid-flag prediction in the fixed reference cohort fails complete-coverage acceptance. Do not filter events or replace outputs.',
        reference_uncertainty='RICH beta is used as the specified high-resolution evaluation reference, not proven exact truth; bootstrap excludes common reference-systematic and upstream reconstruction uncertainties',
        width_inference='Two-primary multiplicity adjustment only for the primary beta-W68 endpoint; all remaining metrics, branches and bins are diagnostics, not an unrestricted multiple-testing discovery claim',
        no_training_on_new_events=True, new_RICH_values_read=False))
    print('PROSPECTIVE_SCORE_PROTOCOL', len(NOMINAL), 'nominal transported branches', len(pairs), 'comparisons', flush=True)


def check():
    gate()
    p = json.loads((SCORE / 'protocol.json').read_text())
    for path, key in [(Path(__file__), 'source_sha256'), (ROOT / 'analysis/evaluate_additional_mc_20261013.py', 'metric_source_sha256'),
                      (ROOT / 'analysis/reference_prospective_flight.py', 'reference_implementation_sha256'),
                      (OUT / 'prediction_protocol.json', 'prediction_protocol_sha256')]:
        if digest(path) != p[key]:
            raise ValueError('Prospective score input changed: ' + key)
    return p


def prediction_coverage(predictions, name, cohort):
    value = np.asarray(predictions[name], float)[cohort]
    valid = np.isfinite(value) & (value != 0)
    if name + '_valid' in predictions:
        valid &= np.asarray(predictions[name + '_valid'], bool)[cohort]
    return dict(complete=bool(valid.all()), invalid_events=int((~valid).sum()), events=len(value))


def score():
    p = check()
    if (SCORE / 'evaluation.json').exists():
        raise FileExistsError('Prospective evaluation already saved')
    audit = json.loads((REF / 'audit.json').read_text())
    if not audit['complete'] or digest(REF / 'reference.npz') != audit['reference_sha256']:
        raise ValueError('Prospective native RICH audit incomplete or changed')
    a = load(OUT / 'measured.npz')
    pred = load(OUT / 'predictions_before_reference.npz')
    ref = load(REF / 'reference.npz')
    for key in ['run', 'event', 'source', 'entry']:
        if not np.array_equal(a[key], pred[key]) or not np.array_equal(a[key], ref[key]):
            raise ValueError('Prospective score identity mismatch: ' + key)
    cohort = ref['valid'].astype(bool)
    n = int(cohort.sum())
    if n < 20:
        raise ValueError('Fewer than twenty valid prospective reference events; preserve the cohort and report insufficient evaluation precision')
    reference_beta = ref['beta'][cohort].astype(float)
    if not np.isfinite(reference_beta).all() or np.any(reference_beta <= 0):
        raise ValueError('Invalid accepted reference beta')
    runs = a['run'][cohort]
    groups = [np.flatnonzero(runs == r) for r in np.unique(runs)]
    if len(groups) < 2:
        raise ValueError('At least two evaluation runs required for run resampling')
    names = sorted(k for k in pred if '_RICH_' in k and k.endswith('_u')) + CONTROLS
    coverage = {name: prediction_coverage(pred, name, cohort) for name in names}
    accepted = [name for name in names if coverage[name]['complete']]
    values = np.column_stack([pred[name][cohort] for name in accepted])
    errors = dict(beta=1 / values - reference_beta[:, None], u=values - 1 / reference_beta[:, None])
    point = {space: matrix_metrics(error) for space, error in errors.items()}
    metrics = {space: {name: dict(zip(METRICS, point[space][j])) for j, name in enumerate(accepted)} for space in errors}
    sensitivities = {}
    for space in errors:
        rows = {}
        for name in accepted:
            if name in CONTROLS:
                continue
            family = next(f for f in FAMILIES if name.startswith(f + '_'))
            tail = name[len(family) + 1:]
            mode_start = tail.index('_fit_RICH_') if '_fit_RICH_' in tail else tail.index('_unit_RICH_')
            branch = tail[:mode_start]
            baseline = family + '_ols' + tail[mode_start:]
            if branch == 'ols' or baseline not in accepted:
                continue
            change = matrix_difference(point[space][accepted.index(baseline)], point[space][accepted.index(name)])
            rows[name] = dict(baseline=baseline, point=dict(zip(DELTA, change)))
        sensitivities[space] = rows
    nominal = [name for name in p['nominal_outputs'] if coverage[name]['complete']]
    nominal_index = [accepted.index(name) for name in nominal]
    nominal_error = {space: error[:, nominal_index] for space, error in errors.items()}
    bootstrap = p['bootstrap']
    count = int(bootstrap['replicates'])
    boot = {space: np.empty((count, len(nominal), len(METRICS))) for space in errors}
    rng = np.random.RandomState(bootstrap['seed'])
    for rep in range(count):
        ix = np.concatenate([groups[j] for j in rng.randint(0, len(groups), len(groups))])
        for space in errors:
            boot[space][rep] = matrix_metrics(nominal_error[space][ix])
        if (rep + 1) % 512 == 0:
            print('PROSPECTIVE_EVALUATION_BOOTSTRAP', rep + 1, '/', count, flush=True)
    intervals = {}
    for space in errors:
        rows = {}
        for new, old in p['pairs']:
            label = new + '_versus_' + old
            if new not in nominal or old not in nominal:
                rows[label] = dict(complete=False, reason='Candidate or baseline failed complete cohort coverage')
                continue
            change = matrix_difference(point[space][accepted.index(old)], point[space][accepted.index(new)])
            samples = matrix_difference(boot[space][:, nominal.index(old)], boot[space][:, nominal.index(new)])
            rows[label] = dict(complete=True, point=dict(zip(DELTA, change)),
                paired_run_interval={key: np.percentile(samples[:, j], bootstrap['quantiles']) for j, key in enumerate(DELTA)},
                two_primary_adjusted_interval={key: np.percentile(samples[:, j], bootstrap['two_primary_adjusted_quantiles']) for j, key in enumerate(DELTA)})
        intervals[space] = rows
    nominal_point = {space: point[space][nominal_index] for space in errors}
    coordinates = dict(measured_R=a['R'][cohort], reference_beta=reference_beta,
        year=np.array([datetime.utcfromtimestamp(int(run)).year for run in runs]))
    bins = []
    for key, edges in p['bins'].items():
        for low, high in zip(edges[:-1], edges[1:]):
            mask = (coordinates[key] >= low) & (coordinates[key] < high)
            if mask.sum() < p['minimum_bin_events']:
                bins.append(dict(coordinate=key, range=[low, high], events=int(mask.sum()), enough_events=False))
                continue
            row = dict(coordinate=key, range=[low, high], events=int(mask.sum()), runs=len(np.unique(runs[mask])), enough_events=True, spaces={})
            for space in errors:
                met = matrix_metrics(nominal_error[space][mask])
                comparisons = {}
                for new, old in p['primary_comparisons']:
                    if new in nominal and old in nominal:
                        comparisons[new] = dict(baseline=old, point=dict(zip(DELTA, matrix_difference(met[nominal.index(old)], met[nominal.index(new)]))))
                row['spaces'][space] = dict(metrics={name: dict(zip(METRICS, met[j])) for j, name in enumerate(nominal)}, primary_comparisons=comparisons)
            bins.append(row)
    primary = []
    for new, old in p['primary_comparisons']:
        label = new + '_versus_' + old
        row = intervals['beta'].get(label)
        supported = bool(row and row['complete'] and row['two_primary_adjusted_interval']['w68_improvement'][0] > 0)
        primary.append(dict(candidate=new, baseline=old, complete_coverage=coverage[new]['complete'] and coverage[old]['complete'],
            multiplicity_adjusted_beta_W68_improvement_supported=supported,
            result=row, scope='This frozen candidate and fixed new NaF reference cohort; not a global optimum or proof that RICH is noiseless'))
    np.savez_compressed(SCORE / 'bootstrap_metrics.npz', **boot)
    result = dict(complete=True, created_utc=stamp(), measured_events=len(a['run']), reference_events=n, runs=len(groups),
        coverage=coverage, metrics=metrics, material_sensitivity=sensitivities, nominal_outputs=nominal,
        comparisons=intervals, primary=primary, bins=bins,
        beta_above_one_retained={name: int(((pred[name][cohort] > 0) & (pred[name][cohort] < 1)).sum()) for name in accepted},
        negative_inverse_speed_retained={name: int((pred[name][cohort] < 0).sum()) for name in accepted},
        reference_error_percentiles=np.percentile(ref['error'][cohort], [0, 5, 50, 95, 100]),
        output_freeze_sha256=digest(OUT / 'output_freeze.json'), reference_sha256=audit['reference_sha256'],
        protocol_sha256=digest(SCORE / 'protocol.json'), bootstrap_metrics_sha256=digest(SCORE / 'bootstrap_metrics.npz'),
        no_output_based_event_removal=True, refit_on_new_events=False, MC_coefficients_used=False)
    write(SCORE / 'evaluation.json', result)
    for row in primary:
        print('PROSPECTIVE_PRIMARY', row['candidate'], json.dumps(row, default=lambda v: v.tolist() if hasattr(v, 'tolist') else v), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['plan', 'score'])
    globals()[parser.parse_args().stage]()
