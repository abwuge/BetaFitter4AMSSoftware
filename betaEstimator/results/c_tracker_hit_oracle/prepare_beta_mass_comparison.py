#!/usr/bin/env python3
"""Join frozen beta and rigidity predictions on the independent tracker test set."""
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'c_joint_beta_rigidity_closure'))
from plot_beta_residual_trends import fit_core
from summarize_strict import load_stage, w68, auc


def save_table(path, rows):
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    base_path = HERE.parent / 'c_mass_toy_oracle' / 'predictions.npz'
    beta_path = HERE.parent / 'c_joint_beta_rigidity_closure' / 'expanded_validation' / 'predictions.npz'
    with np.load(str(base_path)) as cache:
        indices = np.flatnonzero(cache['role'] == b'test')
        data = {k: cache[k][indices] for k in
                ['run', 'event', 'label', 'source', 'u', 'R', 'u_true', 'R_true']}
    assert len(indices) == 20775
    with np.load(str(beta_path)) as cache:
        identities = list(zip(cache['label'], cache['run'], cache['event']))
        mapping = {tuple(int(v) for v in key): i for i, key in enumerate(identities)}
        assert len(mapping) == len(identities)
        matched = np.array([mapping[(int(l), int(r), int(e))] for l, r, e in
                            zip(data['label'], data['run'], data['event'])])
        for key, other in [('u', 'u_raw'), ('R', 'R_raw'),
                           ('u_true', 'u_true'), ('R_true', 'R_true')]:
            assert np.array_equal(data[key], cache[other][matched]), key
        assert cache['truth_valid'][matched].all()
        u_new = cache['u_beta_only_frozen'][matched]
    prediction = load_stage('strict_geometry')
    prediction_indices = prediction['event_index'].astype(np.int64)
    assert len(np.unique(prediction_indices)) == len(indices)
    assert np.array_equal(np.sort(prediction_indices), indices)
    mapping = dict(zip(prediction_indices, prediction['BDTG']))
    correction = np.array([mapping[int(i)] for i in indices], dtype=float)
    assert np.all(np.isfinite(correction)) and np.all(np.isfinite(u_new))
    r_new = data['R'] * np.exp(-correction)
    truth = 1 / data['u_true']
    masses = np.array([11.17793, 12.11221])
    beta_rows, mass_rows, metrics = [], [], []
    edges = np.linspace(.68, .98, 16)
    mass_values = []
    for method, u, rigidity in [(0, data['u'], data['R']), (1, u_new, r_new)]:
        residual = 1000 * (1 / u - truth)
        for low, high in zip(edges[:-1], edges[1:]):
            mask = (truth >= low) & (truth < high)
            if mask.sum() < 80:
                continue
            beta_rows.append(dict(method=method, beta_center=(low + high) / 2,
                                  n=int(mask.sum()), **fit_core(residual[mask])))
        mass2 = (6 * rigidity)**2 * (u**2 - 1)
        mass = np.full(len(u), -1., dtype=float)
        positive = mass2 > 0
        mass[positive] = np.sqrt(mass2[positive])
        mass_values.append(mass)
        for isotope in (0, 1):
            mask = data['label'] == isotope
            fit_mask = mask & positive
            mass_rows.append(dict(isotope=isotope, method=method, n=int(mask.sum()),
                                  n_positive=int(fit_mask.sum()), **fit_core(mass[fit_mask])))
            metrics.append(dict(isotope=isotope, method=method, n=int(mask.sum()),
                                beta_residual_w68=w68(residual[mask]) / 1000,
                                mass_w68=w68(mass[fit_mask]),
                                mass2_error_w68=w68((mass2 - masses[data['label']]**2)[mask]),
                                nonpositive_mass2=int(np.sum(mask & ~positive))))
        metrics.append(dict(isotope=-1, method=method, n=len(u),
                            beta_residual_w68=w68(residual) / 1000,
                            mass_w68=w68(mass[positive]),
                            mass2_error_w68=w68(mass2 - masses[data['label']]**2),
                            nonpositive_mass2=int(np.sum(~positive)),
                            mass2_auc=auc(mass2, data['label'])))
    save_table(HERE / 'beta_mass_beta_fits.tsv', beta_rows)
    save_table(HERE / 'beta_mass_peak_fits.tsv', mass_rows)
    # The text input keeps ROOT plotting independent of Python library versions.
    save_table(HERE / 'beta_mass_events.tsv', [dict(
        isotope=int(data['label'][i]), mass_raw=mass_values[0][i], mass_corrected=mass_values[1][i])
        for i in range(len(indices))])
    (HERE / 'beta_mass_comparison.json').write_text(json.dumps(dict(
        n=len(indices), C12=int(sum(data['label'] == 0)), C13=int(sum(data['label'] == 1)),
        beta_new='1/u_beta_only_frozen; frozen beta-only',
        rigidity_new='R_raw * exp(-strict_geometry BDTG prediction); MC-truth oracle',
        mass='6 R sqrt(u^2-1) [GeV/c^2]',
        beta_residual='beta_reco-beta_true; fit table units 1e-3',
        beta_edges=edges.tolist(), minimum_beta_bin_events=80,
        beta_events_in_range=int(sum(r['n'] for r in beta_rows if r['method'] == 0)),
        truth='Tracker L5 mevmom1[10]; reconstructed Kalman rigidity at z=0',
        selection='Original reconstructed rigidity 2.5-6 GV; whole-source-file tracker test split',
        fit='Conditional unbinned Gaussian MLE, median +/-2.5 W68; observed-Hessian 1sigma errors',
        histogram='Raw counts, 0.1 GeV/c^2 bins; nonpositive mass2 stored in underflow',
        source_hashes={str(p.relative_to(HERE.parent)): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in [base_path, beta_path]}, metrics=metrics), indent=2) + '\n')
    print(json.dumps(metrics, indent=2))
    print('Matched %d test events; %d beta events in plotted range' %
          (len(indices), sum(r['n'] for r in beta_rows if r['method'] == 0)))


if __name__ == '__main__':
    main()
