#!/usr/bin/env python3
"""Compare frozen models on the same flight RICH reference events."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import train_flight_mass_beta as tr
import evaluate_rich_transport as ev

COLORS = ['#5790fc', '#f89c20']  # ROOT Petroff P6 blue and yellow.


def histogram(dest, name, first, second, labels, n):
    fig, ax = plt.subplots(figsize=(8.8, 6.2))
    edges = np.arange(-40, 40.1, 2.)
    for values, label, color in zip((first, second), labels, COLORS):
        values = values * 1000
        ax.hist(values, bins=edges, histtype='step', linewidth=2., color=color,
                label='%s\n$W_{68}$ = %.2f; median = %+.2f' % (label, ev.width(values), np.median(values)))
    ax.set_xlabel(r'$(\widehat\beta_{\rm RICH}-\beta_{\rm RICH}^{\rm meas})\times10^3$', fontsize=16)
    ax.set_ylabel('Events / 0.002', fontsize=16)
    ax.set_xlim(-40, 40)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.30)
    ax.tick_params(labelsize=13)
    ax.legend(frameon=False, loc='upper center', ncol=2, fontsize=11.5)
    ax.text(.97, .71, 'Flight holdout\nN = %d' % n, transform=ax.transAxes, ha='right', fontsize=12)
    fig.tight_layout()
    fig.savefig(str(dest / (name + '.pdf')))
    fig.savefig(str(dest / (name + '.png')), dpi=170)
    plt.close(fig)


def main():
    dest = tr.DEST / 'figures'
    dest.mkdir(exist_ok=True)
    with np.load(str(tr.DEST / 'holdout_predictions.npz')) as f:
        data = {k: f[k] for k in f.files}
    mask = data['reference_mask']
    ref = data['rich_beta'][mask]
    raw = data['raw_beta_rich'][mask] - ref
    flight = data['nominal_beta_rich'][mask] - ref
    mc = data['old_MC_beta_rich'][mask] - ref
    histogram(dest, 'flight_vs_linear', raw, flight, ['Linear TOF', 'Flight mass fit'], len(ref))
    histogram(dest, 'flight_vs_MC', mc, flight, ['MC weak model', 'Flight mass fit'], len(ref))
    result = json.loads((tr.DEST / 'results.json').read_text())
    fig, ax = plt.subplots(figsize=(8.8, 5.4))
    rows = [('MC weak model', result['prior_MC_model']),
            ('Flight: nominal', result['models']['nominal']),
            ('Flight: fraction -0.05', result['models']['fraction_minus_0.05']),
            ('Flight: fraction +0.05', result['models']['fraction_plus_0.05'])]
    for i, (label, values) in enumerate(rows):
        ci = values['bootstrap']['beta_gain_pct_q025_q50_q975']
        color = COLORS[0] if i == 0 else COLORS[1]
        ax.plot([ci[0], ci[2]], [i, i], color=color, linewidth=2.)
        ax.plot(values['gain_beta_pct'], i, 'o', color=color, markersize=7)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([row[0] for row in rows], fontsize=12)
    ax.invert_yaxis()
    ax.set_xlabel(r'Reduction of $\beta$ residual $W_{68}$ vs linear TOF (%)', fontsize=14)
    ax.tick_params(axis='x', labelsize=12)
    ax.text(.98, .94, 'Flight holdout: N = %d\nPaired source bootstrap: 95%% intervals' % len(ref),
            transform=ax.transAxes, ha='right', va='top', fontsize=11)
    fig.tight_layout()
    fig.savefig(str(dest / 'width_gain.pdf'))
    fig.savefig(str(dest / 'width_gain.png'), dpi=170)
    plt.close(fig)
    print('FIGURES COMPLETE', flush=True)


if __name__ == '__main__':
    main()
