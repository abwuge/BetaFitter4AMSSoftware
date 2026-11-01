#!/usr/bin/env python3
"""Plot the saved anchored-energy probe on the expanded validation sample."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'expanded_validation' / 'energy_comparison'
METHODS = [
    ('raw', 'Raw linear TOF', '#5790fc'),
    ('joint_anchored', 'Anchored joint', '#f89c20'),
    ('beta_only_frozen', 'Beta-only frozen', '#e42536'),
    ('anchored_energy', 'Anchored energy', '#964a8b'),
]


def main():
    result = json.loads((OUT / 'probe_results.json').read_text())
    data = result['expanded']
    rows = data['beta_trend']
    grouped = {m: sorted([r for r in rows if r['method'] == m],
                         key=lambda r: r['beta_center']) for m, _, _ in METHODS}
    raw_sigma = np.array([r['sigma'] for r in grouped['raw']])
    centers = np.array([r['beta_center'] for r in grouped['raw']])
    for method, _, _ in METHODS:
        assert len(grouped[method]) == 15
        assert np.allclose(centers, [r['beta_center'] for r in grouped[method]])
        assert [r['n'] for r in grouped[method]] == [r['n'] for r in grouped['raw']]
        for r, sigma in zip(grouped[method], raw_sigma):
            r['sigma_gain_pct'] = 100 * (1 - r['sigma'] / sigma)
    with (OUT / 'beta_trends.tsv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    plt.rcParams.update({'font.size': 13, 'axes.labelsize': 16,
                         'axes.titlesize': 16, 'legend.fontsize': 12,
                         'xtick.direction': 'in', 'ytick.direction': 'in',
                         'xtick.top': True, 'ytick.right': True, 'pdf.fonttype': 42})

    def save(fig, name, pdf):
        fig.tight_layout()
        pdf.savefig(fig)
        fig.savefig(OUT / (name + '.png'), dpi=180)
        fig.savefig(OUT / (name + '.pdf'))
        plt.close(fig)

    with PdfPages(OUT / 'beta_energy_comparison.pdf') as pdf:
        for metric, ylabel, title, name in [
            ('sigma', r'$\sigma(\beta-\beta_{\rm truth})\;[10^{-3}]$',
             'C12 + C13: Gaussian core resolution', 'beta_sigma'),
            ('sigma_gain_pct', r'$1-\sigma_{\rm method}/\sigma_{\rm raw}$ [%]',
             'C12 + C13: resolution improvement', 'beta_sigma_gain'),
            ('mean', r'$\mu(\beta-\beta_{\rm truth})\;[10^{-3}]$',
             'C12 + C13: Gaussian core bias', 'beta_mean'),
        ]:
            fig, ax = plt.subplots(figsize=(11.5, 6.5))
            for method, label, color in METHODS:
                if metric == 'sigma_gain_pct' and method == 'raw':
                    continue
                error = 'sigma_gain_error' if metric == 'sigma_gain_pct' else metric + '_bootstrap_error'
                ax.errorbar(centers, [r[metric] for r in grouped[method]],
                            yerr=[r[error] for r in grouped[method]],
                            fmt='o-', color=color, label=label, lw=1.8, ms=4.5,
                            elinewidth=1.2, capsize=2.5)
            ax.set(xlabel=r'$\beta_{\rm truth}$', ylabel=ylabel,
                   title=title, xlim=(.68, .98))
            ax.margins(y=.25)
            if metric == 'sigma':
                ax.set_ylim(0, max(r['sigma'] for r in rows) * 1.38)
            if metric == 'sigma_gain_pct':
                low = min(r[metric] - r['sigma_gain_error'] for r in rows)
                high = max(r[metric] + r['sigma_gain_error'] for r in rows)
                ax.set_ylim(min(-1, low - .5), high + 3)
            if metric == 'mean':
                ax.set_ylim(min(r['mean'] - r['mean_bootstrap_error'] for r in rows) - .15,
                            max(r['mean'] + r['mean_bootstrap_error'] for r in rows) + .55)
            legend = ax.legend(frameon=False, loc='upper left', ncol=2)
            shown = METHODS[1:] if metric == 'sigma_gain_pct' else METHODS
            for text, (_, _, color) in zip(legend.get_texts(), shown):
                text.set_color(color)
            save(fig, name, pdf)

        fig, ax = plt.subplots(figsize=(11.5, 6.5))
        for i, (method, label, color) in enumerate(METHODS):
            method_result = data['summary'][method]
            gain = method_result['du_gain_vs_raw_pct']
            ax.plot(gain, i, 'o', color=color, ms=9)
            if method != 'raw':
                interval = result['bootstrap']['summary'][method]
                ax.hlines(i, interval['gain_lo95'], interval['gain_hi95'], color=color, lw=3)
                ax.text(gain, i - .23, '95%% CI: %.2f–%.2f%%' %
                        (interval['gain_lo95'], interval['gain_hi95']), ha='center',
                        color=color, fontsize=12)
            ax.text(7.05, i, '%.2f%%    %.5f' % (gain, method_result['w68_du'] * 1000),
                    color=color, va='center', fontsize=13)
        ax.text(7.05, -.62, r'Gain      $W_{68}\;[10^{-3}]$', fontsize=13)
        ax.set_yticks(range(4))
        ax.set_yticklabels([label for _, label, _ in METHODS])
        for label, (_, _, color) in zip(ax.get_yticklabels(), METHODS):
            label.set_color(color)
        ax.set(xlim=(-.4, 10.7), ylim=(3.65, -.9),
               xlabel=r'$1-W_{68,\,\rm method}/W_{68,\,\rm raw}$ [%]',
               title=r'C12 + C13: $W_{68}(1/\beta-1/\beta_{\rm truth})$')
        ax.set_xticks([0, 2, 4, 6])
        save(fig, 'inverse_beta_summary', pdf)
    print('Saved four figures and beta_energy_comparison.pdf in', OUT)


if __name__ == '__main__':
    main()
