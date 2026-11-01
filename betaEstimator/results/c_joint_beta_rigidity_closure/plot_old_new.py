#!/usr/bin/env python3
"""Compare absolute beta errors and reconstructed masses using frozen models."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from bootstrap_energy_comparison import fit_core
from plot_beta_residual_trends import EDGES

HERE = Path(__file__).resolve().parent / 'expanded_validation'
OUT = HERE / 'energy_comparison' / 'old_new'
METHODS = [('old', 'Old: linear TOF', '#5790fc', '--'),
           ('new', 'New: Anchored energy', '#f89c20', '-')]
MASS_EDGES = np.linspace(4, 20, 81)
REPLICATES = 300


def moments(values, bins):
    selected = (bins >= 0) & (bins < len(EDGES) - 1)
    b, x = bins[selected], values[selected]
    count = np.bincount(b, minlength=15)
    total = np.bincount(b, weights=x, minlength=15)
    square = np.bincount(b, weights=x*x, minlength=15)
    assert np.all(count > 1)
    mean = total / count
    variance = (square - total*total/count) / (count - 1)
    return count, mean, np.sqrt(np.maximum(variance, 0))


def main():
    OUT.mkdir(exist_ok=True)
    with np.load(HERE / 'predictions.npz') as f:
        data = dict(f)
    with np.load(HERE / 'energy_comparison' / 'energy_predictions.npz') as f:
        for key in ['run', 'event', 'label', 'source']:
            assert np.array_equal(data[key], f[key]), key
        new_u = f['u_anchored_energy']
    good = data['truth_valid'] & np.isfinite(data['u_true']) & (data['u_true'] > 0)
    u = np.array([data['u_raw'][good], new_u[good]])
    rigidity = data['R_raw'][good]
    assert np.array_equal(data['R_raw'], data['R'])
    assert np.all(np.isfinite(u)) and np.all(u > 0)
    truth = 1 / data['u_true'][good]
    absolute = 1000 * np.abs(1/u - truth)
    labels, sources = data['label'][good], data['source'][good]
    beta_bins = np.searchsorted(EDGES, truth, side='right') - 1
    mass2 = (6*rigidity)**2 * (u*u - 1)
    physical = np.all(np.isfinite(mass2) & (mass2 >= 0), axis=0)
    mass = np.sqrt(np.maximum(mass2, 0))
    groups = {label: [np.flatnonzero((labels == label) & (sources == source))
                      for source in np.unique(sources[labels == label])]
              for label in (0, 1)}
    full_moments = [moments(x, beta_bins) for x in absolute]
    full_fits = [[fit_core(mass[m, physical & (labels == label)])
                  for label in (0, 1)] for m in (0, 1)]
    moment_draws = np.empty((REPLICATES, 2, 15, 2))
    mass_draws = np.empty((REPLICATES, 2, 2, 2))
    rng = np.random.RandomState(20261007)
    for rep in range(REPLICATES):
        idx = np.concatenate([groups[label][j] for label in (0, 1)
                              for j in rng.randint(len(groups[label]), size=len(groups[label]))])
        for m in (0, 1):
            _, average, spread = moments(absolute[m, idx], beta_bins[idx])
            moment_draws[rep, m, :, 0] = average
            moment_draws[rep, m, :, 1] = spread
            for label in (0, 1):
                selected = idx[physical[idx] & (labels[idx] == label)]
                fit = fit_core(mass[m, selected])
                mass_draws[rep, m, label] = [fit['mean'], fit['sigma']]
        if (rep + 1) % 100 == 0:
            print(f'Paired source-file bootstrap {rep + 1}/{REPLICATES}', flush=True)
    assert np.all(np.isfinite(moment_draws)) and np.all(np.isfinite(mass_draws))
    moment_errors = np.std(moment_draws, axis=0, ddof=1)
    mass_errors = np.std(mass_draws, axis=0, ddof=1)
    rows = []
    centers = (EDGES[:-1] + EDGES[1:])/2
    for m, (method, _, _, _) in enumerate(METHODS):
        count, average, spread = full_moments[m]
        for b, center in enumerate(centers):
            rows.append({'method': method, 'beta_low': EDGES[b], 'beta_high': EDGES[b+1],
                         'beta_center': center, 'n': int(count[b]),
                         'mean_absolute_residual': float(average[b]),
                         'mean_bootstrap_error': float(moment_errors[m, b, 0]),
                         'std_absolute_residual': float(spread[b]),
                         'std_bootstrap_error': float(moment_errors[m, b, 1])})
    with (OUT / 'beta_absolute_trends.tsv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    hist_rows, mass_summary = [], []
    for label in (0, 1, -1):
        selected = physical & ((labels == label) if label >= 0 else True)
        sample = f'C{12+label}' if label >= 0 else 'combined'
        for m, (method, _, _, _) in enumerate(METHODS):
            values = mass[m, selected]
            q16, median, q84 = np.percentile(values, [16, 50, 84])
            counts, _ = np.histogram(values, MASS_EDGES)
            summary = {'sample': sample, 'method': method, 'n': len(values),
                       'median': float(median), 'w68': float((q84-q16)/2),
                       'below_plot': int(np.sum(values < MASS_EDGES[0])),
                       'above_plot': int(np.sum(values > MASS_EDGES[-1]))}
            if label >= 0:
                summary.update(full_fits[m][label])
                summary['mean_error'], summary['sigma_error'] = mass_errors[m, label].tolist()
            mass_summary.append(summary)
            for b, count in enumerate(counts):
                hist_rows.append({'sample': sample, 'method': method,
                                  'mass_low': MASS_EDGES[b], 'mass_high': MASS_EDGES[b+1],
                                  'count': int(count)})
    with (OUT / 'mass_histograms.tsv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(hist_rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(hist_rows)
    plt.rcParams.update({'font.size': 13, 'axes.labelsize': 16, 'axes.titlesize': 16,
                         'legend.fontsize': 12, 'xtick.direction': 'in',
                         'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
                         'pdf.fonttype': 42})
    pages = []

    def finish(fig, ax, name, pdf):
        legend = ax.legend(frameon=False, loc='upper left')
        for text, (_, _, color, _) in zip(legend.get_texts(), METHODS):
            text.set_color(color)
        fig.tight_layout()
        pdf.savefig(fig)
        fig.savefig(OUT / (name + '.pdf'))
        fig.savefig(OUT / (name + '.png'), dpi=180)
        plt.close(fig)
        pages.append(name)

    with PdfPages(OUT / 'old_new_comparison.pdf') as pdf:
        for metric, ylabel, title, name in [
            (0, r'$\langle |\beta_{\rm reco}-\beta_{\rm truth}|\rangle\;[10^{-3}]$',
             'C12 + C13: mean absolute beta residual', 'beta_absolute_mean'),
            (1, r'$\mathrm{SD}(|\beta_{\rm reco}-\beta_{\rm truth}|)\;[10^{-3}]$',
             'C12 + C13: absolute beta residual spread', 'beta_absolute_std')]:
            fig, ax = plt.subplots(figsize=(10, 6.5))
            upper = 0
            for m, (_, label, color, style) in enumerate(METHODS):
                values = full_moments[m][metric+1]
                errors = moment_errors[m, :, metric]
                ax.errorbar(centers, values, yerr=errors, fmt='o', linestyle=style,
                            color=color, lw=2, ms=4.5, capsize=3, label=label)
                upper = max(upper, np.max(values+errors))
            ax.set(xlabel=r'$\beta_{\rm truth}$', ylabel=ylabel, title=title,
                   xlim=(.68, .98), ylim=(0, upper*1.28))
            finish(fig, ax, name, pdf)
        for isotope in (0, 1, -1):
            sample = f'C{12+isotope}' if isotope >= 0 else 'combined'
            selected = physical & ((labels == isotope) if isotope >= 0 else True)
            fig, ax = plt.subplots(figsize=(10, 6.5))
            ymax = 0
            for m, (method, label, color, style) in enumerate(METHODS):
                counts, _ = np.histogram(mass[m, selected], MASS_EDGES)
                ymax = max(ymax, counts.max())
                summary = next(r for r in mass_summary if r['sample'] == sample and r['method'] == method)
                if isotope >= 0:
                    fit_label = (label + '\n' +
                                 r'$\mu=%.3f\pm%.3f,\;\sigma=%.3f\pm%.3f$ GeV/$c^2$' %
                                 (summary['mean'], summary['mean_error'], summary['sigma'], summary['sigma_error']))
                else:
                    fit_label = label + '\n' + r'median $=%.3f$, $W_{68}=%.3f$ GeV/$c^2$' % (summary['median'], summary['w68'])
                ax.hist(mass[m, selected], bins=MASS_EDGES, histtype='step',
                        color=color, linestyle=style, lw=2.2, label=fit_label)
            title = f'{sample}: reconstructed mass' if isotope >= 0 else 'C12 + C13: observed combined mass spectrum'
            ax.set(xlabel=r'Reconstructed mass [GeV/$c^2$]',
                   ylabel=r'Events / (0.2 GeV/$c^2$)', title=title,
                   xlim=(4, 20), ylim=(0, ymax*1.55))
            ax.text(.97, .95, f'N = {selected.sum():,}\n2.5 < R < 6 GV',
                    transform=ax.transAxes, ha='right', va='top', fontsize=12)
            finish(fig, ax, 'mass_' + sample, pdf)
    metadata = {'old': 'raw linear TOF', 'new': 'frozen Anchored energy',
                'residual': 'abs(beta_reco - beta_truth)', 'residual_units': '1e-3',
                'beta_statistics': 'untrimmed arithmetic mean and sample standard deviation; no Gaussian fit',
                'truth_valid_events': int(good.sum()),
                'beta_events_in_range': int(full_moments[0][0].sum()),
                'mass_formula': 'm = 6 R_raw sqrt(u_method^2 - 1), GeV/c^2',
                'mass_common_events': int(physical.sum()),
                'negative_mass2_old': int(np.sum(mass2[0] < 0)),
                'negative_mass2_new': int(np.sum(mass2[1] < 0)),
                'mass_excluded_common': int(np.sum(~physical)),
                'mass_selection': 'same truth-valid sample; require nonnegative mass2 for both old and new',
                'mass_fit': 'per-isotope conditional Gaussian core in median +/- 2.5 W68; combined spectrum unfitted',
                'bootstrap': {'replicates': REPLICATES, 'seed': 20261007,
                              'unit': 'source file, stratified by isotope, paired between old/new',
                              'error_bars': 'one bootstrap standard deviation, fixed models'},
                'mass_summary': mass_summary, 'pages': pages}
    (OUT / 'summary.json').write_text(json.dumps(metadata, indent=2) + '\n')
    np.savez_compressed(OUT / 'bootstrap.npz', beta_moments=moment_draws, mass_fits=mass_draws)
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == '__main__':
    main()
