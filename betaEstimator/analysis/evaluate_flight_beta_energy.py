#!/usr/bin/env python3
"""Energy-dependent diagnostics of the unchanged flight beta predictions."""
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import train_flight_mass_beta as tr
import evaluate_rich_transport as ev
from evaluate_flight_mass_beta import check_freeze

MASS_PER_NUCLEON = 11.17793 / 12.
BETA_EDGES = np.array([.78, .85, .90, .95, 1.02])
COLORS = ['#5790fc', '#f89c20']


def kinetic_energy(beta):
    beta = np.asarray(beta, float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return MASS_PER_NUCLEON * (1 / np.sqrt(1 - beta * beta) - 1)


def paired_errors(old, new, reference, source):
    groups = [np.flatnonzero(source == s) for s in np.unique(source)]
    rng = np.random.RandomState(20261007)
    values = []
    for _ in range(400):
        ix = np.concatenate([groups[i] for i in rng.randint(0, len(groups), len(groups))])
        r0, r1 = old[ix] - reference[ix], new[ix] - reference[ix]
        values.append([ev.width(r0), ev.width(r1), np.mean(r0), np.mean(r1),
                       100 * (1 - ev.width(r1) / ev.width(r0))])
    values = np.array(values)
    return {'repeats': 400, 'source_groups': len(groups),
            'old_W68_sd': np.std(values[:, 0], ddof=1),
            'new_W68_sd': np.std(values[:, 1], ddof=1),
            'old_mean_sd': np.std(values[:, 2], ddof=1),
            'new_mean_sd': np.std(values[:, 3], ddof=1),
            'gain_95ci': np.percentile(values[:, 4], [2.5, 97.5])}


def axes():
    fig, ax = plt.subplots(figsize=(8.8, 6.2))
    ax.set_xlim(.5, 2.65)
    ax.set_xlabel(r'$T_{\rm RICH}/A$ (GeV/n), from reference $\beta$', fontsize=15)
    ax.tick_params(labelsize=12)
    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    beta_ticks = np.array([.80, .85, .90, .95, .96])
    top.set_xticks(kinetic_energy(beta_ticks))
    top.set_xticklabels(['%.2f' % b for b in beta_ticks], fontsize=12)
    top.set_xlabel(r'$\beta_{\rm RICH}$', fontsize=15, labelpad=8)
    return fig, ax


def save(fig, name):
    fig.tight_layout()
    folder = tr.DEST / 'figures'
    fig.savefig(str(folder / (name + '.pdf')))
    fig.savefig(str(folder / (name + '.png')), dpi=170)
    plt.close(fig)


def plots(rows):
    x = np.array([r['energy_median_GeV_per_nucleon'] for r in rows])
    for metric, label, filename, error in [
        ('w68', r'$W_{68}(\widehat\beta_{\rm RICH}-\beta_{\rm RICH})\times10^3$', 'energy_width', 'W68'),
        ('mean', r'$\langle\widehat\beta_{\rm RICH}-\beta_{\rm RICH}\rangle\times10^3$', 'energy_bias', 'mean')]:
        fig, ax = axes()
        for method, name, color, marker in [('baseline_beta', 'Linear TOF', COLORS[0], 'o'),
                                           ('new_beta', 'Flight mass fit', COLORS[1], 's')]:
            key = 'old' if method == 'baseline_beta' else 'new'
            ax.errorbar(x, [1000 * r[method][metric] for r in rows],
                        yerr=[1000 * r['errors'][key + '_' + error + '_sd'] for r in rows],
                        color=color, marker=marker, linewidth=1.8, capsize=4, label=name)
        ax.set_ylabel(label, fontsize=14)
        ax.legend(frameon=False, fontsize=12, loc='upper left' if metric == 'w68' else 'lower left')
        ax.text(.97, .97, 'Flight holdout\nErrors: bootstrap SD', transform=ax.transAxes,
                ha='right', va='top', fontsize=11)
        if metric == 'w68':
            ax.set_ylim(0, max(1000 * r[m]['w68'] for r in rows for m in ('baseline_beta', 'new_beta')) * 1.42)
        save(fig, filename)
    fig, ax = axes()
    for r in rows:
        ci = r['errors']['gain_95ci']
        xx = r['energy_median_GeV_per_nucleon']
        ax.plot([xx, xx], ci, color=COLORS[1], linewidth=2.)
        ax.plot(xx, r['gain_beta_pct'], 'o', color=COLORS[1], markersize=7)
        ax.annotate('N = %d' % r['n'], (xx, ci[1]), xytext=(0, 7), textcoords='offset points',
                    ha='center', fontsize=11)
    ax.axhline(0, color=COLORS[0], linewidth=1.2, linestyle='--')
    ax.set_ylabel(r'Reduction of $\beta$ residual $W_{68}$ (%)', fontsize=14)
    ax.text(.97, .05, 'Flight vs linear TOF\nPaired source bootstrap: 95% intervals',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=11)
    ax.set_ylim(-35, 40)
    save(fig, 'energy_gain')


def main():
    check_freeze()
    path = tr.DEST / 'holdout_predictions.npz'
    with np.load(str(path), allow_pickle=False) as f:
        p = {k: f[k] for k in f.files}
    mask = p['reference_mask']
    ref = p['rich_beta']
    energy = kinetic_energy(ref)
    invalid_energy = mask & (~np.isfinite(energy) | (ref >= 1))
    # The original beta bins remain authoritative. Invalid energy references
    # are counted explicitly and never reassigned using a model prediction.
    if invalid_energy.any():
        raise RuntimeError('Selected RICH beta cannot be converted to physical energy')
    original = json.loads((tr.DEST / 'results.json').read_text())['models']['nominal']
    rows = []
    for lo, hi in zip(BETA_EDGES[:-1], BETA_EDGES[1:]):
        sel = mask & (ref >= lo) & (ref < hi)
        old, new = p['raw_beta_rich'][sel], p['nominal_beta_rich'][sel]
        reference = ref[sel]
        result = ev.compare(old, new, reference)
        existing = original['groups']['beta_%.2f_%.2f' % (lo, hi)]
        if result['baseline_beta']['n'] != existing['baseline_beta']['n']:
            raise RuntimeError('Original beta-bin counts changed')
        if not np.isclose(result['gain_beta_pct'], existing['gain_beta_pct'], rtol=0, atol=1e-12):
            raise RuntimeError('Original beta-bin score changed')
        errors = paired_errors(old, new, reference, p['source'][sel])
        expected_ci = np.array(existing['bootstrap']['beta_gain_pct_q025_q50_q975'])[[0, 2]]
        if not np.allclose(errors['gain_95ci'], expected_ci, rtol=0, atol=1e-10):
            raise RuntimeError('Paired source-bootstrap reproduction failed')
        rows.append(dict(result, beta_low=lo, beta_high=hi, n=int(sel.sum()),
                         energy_low_GeV_per_nucleon=float(kinetic_energy(lo)),
                         energy_high_GeV_per_nucleon=float(kinetic_energy(hi)) if hi < 1 else None,
                         energy_median_GeV_per_nucleon=float(np.median(energy[sel])), errors=errors))
    if sum(r['n'] for r in rows) != int(mask.sum()):
        raise RuntimeError('Energy bins do not cover the reference sample')
    ev.dump(tr.DEST / 'energy_results.json', {
        'model_changed': False, 'prediction_sha256': ev.digest(path),
        'script_sha256': ev.digest(Path(__file__)), 'n': int(mask.sum()),
        'energy_definition': 'T_RICH/A = (11.17793/12) * ((1-beta_RICH^2)^(-1/2)-1) GeV/n',
        'mass_per_nucleon_GeV': MASS_PER_NUCLEON,
        'energy_status': 'RICH-plane kinetic energy per nucleon, with fixed C12 m/A approximation; not center-plane truth',
        'binning': 'existing four reference-beta bins, identical for every method',
        'invalid_energy_n': int(invalid_energy.sum()), 'rows': rows})
    with (tr.DEST / 'energy_metrics.tsv').open('w', newline='') as f:
        columns = ['beta_low', 'beta_high', 'energy_low', 'energy_high', 'n', 'raw_W68', 'flight_W68',
                   'gain_pct', 'gain_CI_low', 'gain_CI_high', 'raw_mean', 'flight_mean']
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(columns)
        for r in rows:
            writer.writerow([r['beta_low'], r['beta_high'], r['energy_low_GeV_per_nucleon'],
                             r['energy_high_GeV_per_nucleon'], r['n'], r['baseline_beta']['w68'],
                             r['new_beta']['w68'], r['gain_beta_pct'], *r['errors']['gain_95ci'],
                             r['baseline_beta']['mean'], r['new_beta']['mean']])
    plots(rows)
    print('ENERGY DIAGNOSTICS COMPLETE: unchanged predictions; 710 events in four bins', flush=True)
    for r in rows:
        print(r['beta_low'], r['beta_high'], 'N', r['n'], 'W68 gain', r['gain_beta_pct'],
              '95% CI', r['errors']['gain_95ci'], flush=True)


if __name__ == '__main__':
    main()
