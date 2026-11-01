#!/usr/bin/env python3
"""Bound the effect of observed endpoint coordinates on the frozen primary tests."""
import json
from pathlib import Path
import numpy as np
from reference_prospective_flight_schema import OUT, REF, check
from research_flight_conditional import load, write, digest
from research_conditional_mc import stamp
from evaluate_additional_mc_20261013 import matrix_metrics, matrix_difference, METRICS, DELTA

END = REF / 'endpoint_sensitivity'
SCORE = OUT / 'scoring'


def gain_bounds(old_width, new_width, old_epsilon, new_epsilon):
    old_width = np.asarray(old_width, float); new_width = np.asarray(new_width, float)
    if np.any(old_width <= old_epsilon):
        raise ValueError('The lower endpoint-width denominator is nonpositive')
    lower = 1 - (new_width + new_epsilon) / (old_width - old_epsilon)
    upper = 1 - np.maximum(new_width - new_epsilon, 0) / (old_width + old_epsilon)
    return lower, upper


def main():
    check()
    destination = END / 'score_sensitivity.json'
    if destination.exists():
        raise FileExistsError('Endpoint score sensitivity already saved')
    p = json.loads((END / 'evaluation.json').read_text())
    original = json.loads((SCORE / 'evaluation.json').read_text())
    protocol = json.loads((SCORE / 'protocol.json').read_text())
    if not p['complete'] or not original['complete'] or digest(END / 'predictions.npz') != p['predictions_sha256']:
        raise ValueError('Complete immutable endpoint and primary evaluations required')
    alternate = load(END / 'predictions.npz'); ref = load(REF / 'reference.npz')
    mask = ref['valid'].astype(bool)
    for key in ['run', 'event']:
        if not np.array_equal(alternate[key], ref[key][mask]):
            raise ValueError('Endpoint score cohort identity mismatch')
    names = sorted(key for key in alternate if key.endswith('_u'))
    beta = np.column_stack([1 / alternate[name] for name in names])
    met = matrix_metrics(beta - ref['beta'][mask, None])
    metrics = {name: dict(zip(METRICS, met[j])) for j, name in enumerate(names)}
    if digest(SCORE / 'bootstrap_metrics.npz') != original['bootstrap_metrics_sha256']:
        raise ValueError('Original evaluation-run bootstrap changed')
    with np.load(SCORE / 'bootstrap_metrics.npz') as f:
        boot = f['beta'].copy()
    nominal = original['nominal_outputs']
    rows = []
    for new, old in protocol['primary_comparisons']:
        eo = p['branches'][old]['beta_width_absolute_bound']; en = p['branches'][new]['beta_width_absolute_bound']
        old_metric = original['metrics']['beta'][old]; new_metric = original['metrics']['beta'][new]
        low, high = gain_bounds(old_metric['w68'], new_metric['w68'], eo, en)
        lower, upper = gain_bounds(boot[:, nominal.index(old), METRICS.index('w68')],
                                  boot[:, nominal.index(new), METRICS.index('w68')], eo, en)
        adjusted = protocol['bootstrap']['two_primary_adjusted_quantiles']
        lower_interval = np.percentile(lower, adjusted); upper_interval = np.percentile(upper, adjusted)
        exact = matrix_difference(met[names.index(old)], met[names.index(new)])
        if not (float(low) - 1e-12 <= exact[0] <= float(high) + 1e-12):
            raise ValueError('Observed endpoint result violates the deterministic width bound')
        rows.append(dict(candidate=new, baseline=old, nominal_w68_gain=1-new_metric['w68']/old_metric['w68'],
            observed_endpoint_changes=dict(zip(DELTA, exact)), point_w68_gain_bounds=[float(low), float(high)],
            multiplicity_adjusted_interval_enclosure=[float(lower_interval[0]), float(upper_interval[-1])],
            all_bootstrap_lower_gain_quantiles=lower_interval, all_bootstrap_upper_gain_quantiles=upper_interval,
            adjusted_positive_width_conclusion_invariant=bool(lower_interval[0] > 0),
            deterministic_input_bounds=dict(baseline=eo, candidate=en)))
    write(destination, dict(complete=True, created_utc=stamp(), events=int(mask.sum()), metrics=metrics, primary=rows,
        proof='Coordinate perturbations change every beta residual by at most the recorded epsilon. Empirical quantiles and W68 are 1-Lipschitz in sup norm. Applying interval division to each paired bootstrap width ratio gives a simultaneous deterministic enclosure for that same set of 4096 resamples.',
        interval_interpretation='An enclosure of the original paired-bootstrap endpoint-sensitivity interval; not a fresh independent confidence experiment and not a bound on all material-model errors.',
        original_primary_evaluation_sha256=digest(SCORE / 'evaluation.json'), reference_sha256=p['reference_sha256'],
        endpoint_evaluation_sha256=digest(END / 'evaluation.json'), source_sha256=digest(Path(__file__)),
        coordinate_tolerance_relabelled_as_pass=False, primary_predictions_changed=False, new_hypothesis_selected=False))
    for row in rows:
        print('PROSPECTIVE_ENDPOINT_SCORE', json.dumps(row, default=lambda value: value.tolist() if hasattr(value, 'tolist') else value), flush=True)


if __name__ == '__main__':
    main()
