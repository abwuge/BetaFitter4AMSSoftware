#!/usr/bin/env python3
"""Score frozen flight-only mass fits against the reserved RICH reference."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

import rich_transport as rt
import evaluate_rich_transport as ev
import train_flight_mass_beta as tr


def check_freeze():
    freeze = json.loads((tr.DEST / 'freeze.json').read_text())
    if freeze['code_sha256'] != tr.dependency_hashes():
        raise RuntimeError('Frozen training dependency changed')
    for name, digest in freeze['model_sha256'].items():
        if digest != ev.digest(tr.DEST / (name + '_model.json')):
            raise RuntimeError('Frozen model changed: ' + name)
    if freeze['material_sha256'] != ev.digest(rt.MAP):
        raise RuntimeError('Material map changed')
    return freeze


def prediction_checks():
    """Check serialization and sensitivity to a common timing-slope shift."""
    with np.load(str(tr.CACHE / 'training.npz'), allow_pickle=False) as f:
        d = {k: f[k] for k in f.files}
    checks = {}
    summary = json.loads((tr.DEST / 'training.json').read_text())
    for name in summary['models']:
        model = rt.CenterModel(tr.DEST / (name + '_model.json'))
        matrix = tr.feature_matrix(d, model.state)
        expected = d['u'] - matrix.dot(model.state['parameters'][:matrix.shape[1]])
        actual = model.predict(d)
        difference = float(np.max(np.abs(actual - expected)))
        if difference > 1e-12:
            raise RuntimeError('Serialized model replay failed')
        injection = {}
        for shift in (-.02, -.005, .005, .02):
            shifted = dict(d, u=d['u'] + shift)
            ratio = (model.predict(shifted) - actual) / shift
            injection[str(shift)] = {'median': np.median(ratio),
                                      'q05_q95': np.percentile(ratio, [5, 95]),
                                      'fraction_negative': np.mean(ratio < 0)}
        checks[name] = {'serialization_max_difference': difference,
                        'slope_injection_response': injection,
                        'mean_correction_by_bin': [np.mean((d['u'] - actual)[d['bin'] == b]) for b in range(4)]}
    return checks


def raw_fraction_profiles():
    """Training-spectrum profiles, conditional on the two-component response."""
    with np.load(str(tr.CACHE / 'training.npz'), allow_pickle=False) as f:
        data = {k: f[k] for k in f.files}
    fit = json.loads((tr.DEST / 'training.json').read_text())['raw_spectrum']
    profiles = []
    for b in range(4):
        subset = ev.sub(data, data['bin'] == b)
        density = tr.RawMixture(subset).base
        fractions = np.unique(np.r_[.001, .01, .02, .05, .1, .2, .3, .5, .7, .9, .999, fit['fraction'][b]])
        records = []
        for fraction in fractions:
            density.prior[:] = fraction
            def objective(logsigma):
                params = np.r_[0., np.full(4, logsigma)]
                return np.sum(density.event_nll(params))
            optimized = minimize_scalar(objective, bounds=(np.log(.0015), np.log(.08)), method='bounded')
            if not optimized.success:
                raise RuntimeError('Training fraction profile failed')
            records.append({'fraction': fraction, 'nll_sum': optimized.fun,
                            'sigma_u': np.exp(optimized.x)})
        best = min(r['nll_sum'] for r in records)
        for record in records:
            record['twice_delta_nll'] = 2 * (record['nll_sum'] - best)
        profiles.append({'bin': b, 'n': len(subset['u']), 'scan': records})
    return profiles


def score(old, new, reference, source, mask):
    if mask.sum() < 2:
        return {'n': int(mask.sum())}
    result = ev.compare(old[mask], new[mask], reference[mask])
    result['bootstrap'] = ev.bootstrap(old[mask], new[mask], reference[mask], source[mask])
    result['fraction_absolute_error_reduced'] = np.mean(
        np.abs(new[mask] - reference[mask]) < np.abs(old[mask] - reference[mask]))
    return result


def evaluate():
    freeze = check_freeze()
    verification = prediction_checks()
    profiles = raw_fraction_profiles()
    # Freeze the evaluation implementation and its input identities before
    # opening any of the RICH measurement arrays for this experiment.
    baseline_dir = rt.ROOT / 'results/rich_transport'
    cache = Path('/tmp/betaEstimator-rich-transport/data_evaluation_input.npz')
    inputs = [cache, baseline_dir / 'data_predictions.npz',
              baseline_dir / 'data_reference_v2_masks.npz']
    protocol = {'evaluator_sha256': ev.digest(Path(__file__)),
                'input_sha256': {str(p): ev.digest(p) for p in inputs},
                'nominal_model_sha256': freeze['model_sha256']['nominal'],
                'primary': 'nominal', 'sensitivity_candidates': freeze['sensitivity'],
                'reference': 'existing charge_consistent_reference mask; no beta-agreement cut',
                'metrics': 'W68, mean, median, std, RMSE, p95 absolute; paired source bootstrap 400',
                'groups': 'NaF, aerogel, strict ring, charge 5-7, four RICH-beta ranges',
                'prior_use': 'same holdout inspected previously for the MC-trained model; not a new blind dataset'}
    protocol_path = tr.DEST / 'evaluation_freeze.json'
    if protocol_path.exists() and json.loads(protocol_path.read_text()) != protocol:
        raise RuntimeError('Frozen evaluation protocol changed')
    ev.dump(protocol_path, protocol)
    with np.load(str(cache), allow_pickle=False) as f:
        a = {k: f[k] for k in f.files}
    with np.load(str(baseline_dir / 'data_predictions.npz'), allow_pickle=False) as f:
        baseline = {k: f[k] for k in f.files}
    with np.load(str(baseline_dir / 'data_reference_v2_masks.npz'), allow_pickle=False) as f:
        quality = {k: f[k] for k in f.files}
    for d in (baseline, quality):
        if not np.array_equal(a['run'], d['run']) or not np.array_equal(a['event'], d['event']):
            raise RuntimeError('Holdout event identity mismatch')
    if not np.all(a['outer_holdout']):
        raise RuntimeError('Non-holdout event in RICH evaluation')
    with np.load(str(tr.CACHE / 'training.npz'), allow_pickle=False) as f:
        if set(a['run'].tolist()) & set(f['run'].tolist()):
            raise RuntimeError('Training/holdout run overlap')
    d = rt.center_features(a)
    if np.max(np.abs(d['u'] - baseline['u_raw'])) > 1e-12:
        raise RuntimeError('Raw-beta baseline replay differs')
    centers = {name: rt.CenterModel(tr.DEST / (name + '_model.json')).predict(d)
               for name in freeze['model_sha256']}
    mask = quality['charge_consistent_reference'].astype(bool)
    rows = np.flatnonzero(mask)
    material = rt.MaterialMap()
    geometry = material.build_columns(a['tk_pos'][rows], a['tk_dir'][rows])
    del material
    transport = {name: rt.propagate(u[rows], d['R'][rows], geometry)
                 for name, u in centers.items()}
    # Every candidate, including the previous MC model, uses identical events.
    common_valid = np.ones(len(rows), bool)
    for pred in transport.values():
        common_valid &= pred['valid']
    mask[rows[~common_valid]] = False
    old = baseline['raw_beta_planes'][:, -1]
    reference = a['rich_beta'][:, 0].astype(float)
    q = a['rich_q'][:, 1]
    group_masks = {'NaF': a['rich_NaF'].astype(bool),
                   'aerogel': ~a['rich_NaF'].astype(bool),
                   'strict_ring_quality': baseline['rich_strict'],
                   'charge_5_to_7': (q > 5) & (q < 7)}
    for lo, hi in ((.78, .85), (.85, .90), (.90, .95), (.95, 1.02)):
        group_masks['beta_%.2f_%.2f' % (lo, hi)] = (reference >= lo) & (reference < hi)
    result = {'reference_n': int(mask.sum()), 'existing_reference_n': len(rows),
              'selected_carbon_n': len(reference), 'verification': verification,
              'training_fraction_profiles': profiles,
              'prior_MC_model': score(old, baseline['new_beta_planes'][:, -1], reference, a['source'], mask),
              'models': {}}
    saved = {k: a[k] for k in ('run', 'event', 'source')}
    saved.update(reference_mask=mask, rich_beta=reference, raw_beta_rich=old,
                 old_MC_beta_rich=baseline['new_beta_planes'][:, -1], u_raw=d['u'])
    for name, pred in transport.items():
        new = np.full(len(old), np.nan)
        new[rows] = pred['beta_planes'][:, -1]
        summary = score(old, new, reference, a['source'], mask)
        summary['groups'] = {group: score(old, new, reference, a['source'], mask & g)
                             for group, g in group_masks.items()}
        summary['vs_MC_model'] = score(baseline['new_beta_planes'][:, -1], new, reference, a['source'], mask)
        summary['expected_total_loss_GeV'] = ev.stats(pred['loss_GeV'][common_valid].sum(axis=1))
        result['models'][name] = summary
        saved[name + '_beta_rich'] = new
        saved[name + '_u_center'] = centers[name]
    np.savez_compressed(str(tr.DEST / 'holdout_predictions.npz'), **saved)
    ev.dump(tr.DEST / 'results.json', result)
    ev.dump(tr.DEST / 'verification.json', {'complete': True, 'unique_training_events': True,
            'disjoint_runs': True, 'all_candidates_frozen_before_RICH_access': True,
            'reference_n': int(mask.sum()), 'nominal': verification['nominal']})
    print('FLIGHT EVALUATION COMPLETE', flush=True)
    for name, summary in result['models'].items():
        print(name, 'N', summary['baseline_beta']['n'], 'W68 gain', summary['gain_beta_pct'],
              'CI', summary['bootstrap']['beta_gain_pct_q025_q50_q975'],
              'RMSE', summary['baseline_beta']['rmse'], summary['new_beta']['rmse'], flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.parse_args()
    evaluate()


if __name__ == '__main__':
    main()
