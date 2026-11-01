#!/usr/bin/env python3
"""Compare frozen beta estimates with MC truth for the expanded validation."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

HERE = Path(__file__).resolve().parent / 'expanded_validation'
OUT = HERE / 'figures'


def main():
    with np.load(HERE / 'predictions.npz') as f:
        d = dict(f)
    mask = d['truth_valid'] & np.isfinite(d['u_true']) & (d['u_true'] > 0)
    methods = [
        ('beta_tof', 'u_raw', r'$\beta_{\mathrm{TOF}}$', 'Linear TOF', '#5790fc'),
        ('beta_masspeak', 'u_joint_anchored', r'$\beta_{\mathrm{masspeak}}$', 'Anchored joint', '#f89c20'),
        ('beta_masspeak_beta_only', 'u_beta_only_frozen', r'$\beta_{\mathrm{masspeak}}$', r'Frozen $\beta$-only', '#e42536'),
    ]
    for _, key, _, _, _ in methods:
        assert np.all(np.isfinite(d[key][mask]) & (d[key][mask] > 0))
    x = 1 / d['u_true'][mask]
    assert len(x) == 115091
    # Draw the same deterministic event ordering for every method, without thinning.
    order = np.random.RandomState(13013).permutation(len(x))
    plt.rcParams.update({'font.size': 13, 'axes.labelsize': 17, 'axes.titlesize': 16,
                         'legend.fontsize': 12, 'xtick.direction': 'in', 'ytick.direction': 'in',
                         'xtick.top': True, 'ytick.right': True, 'pdf.fonttype': 42})
    audit = {'events': len(x), 'selection': 'All truth-valid expanded-validation events; no subsampling',
             'main_masspeak': 'u_joint_anchored', 'tof': 'u_raw',
             'core_view': {'xlim': [.65, 1.0], 'ylim': [.65, 1.0]}, 'methods': []}
    with PdfPages(OUT / 'beta_truth_scatter.pdf') as pdf:
        for stem, key, ylabel, method, color in methods:
            y = 1 / d[key][mask]
            outside = int(np.sum((x < .65) | (x > 1) | (y < .65) | (y > 1)))
            q = np.percentile(y-x, [16, 50, 84])
            audit['methods'].append({'method': method, 'prediction_key': key,
                                      'outside_core_view': outside,
                                      'median_beta_residual': q[1], 'w68_beta_residual': (q[2]-q[0])/2})
            for full in (False, True):
                fig, ax = plt.subplots(figsize=(8, 7))
                ax.scatter(x[order], y[order], s=1.5, color=color, alpha=.10,
                           linewidths=0, rasterized=True)
                ax.plot([.4, 1.15], [.4, 1.15], color='0.45', lw=1, ls='--', label=r'$\beta_{\rm reco}=\beta_{\rm truth}$')
                ax.set(xlabel=r'$\beta_{\rm truth}$', ylabel=ylabel,
                       title=method, xlim=(.45, 1.0) if full else (.65, 1.0),
                       ylim=(.45, 1.15) if full else (.65, 1.0))
                ax.set_aspect('equal', adjustable='box')
                ax.text(.035, .96, f'C12 + C13, N = {len(x):,}\n'+r'$2.5 \leq R_{\rm reco} < 6$ GV',
                        transform=ax.transAxes, va='top', fontsize=12)
                ax.legend(loc='lower right', frameon=False)
                fig.tight_layout()
                suffix = '_full' if full else ''
                fig.savefig(OUT / (stem + '_vs_truth' + suffix + '.png'), dpi=190)
                fig.savefig(OUT / (stem + '_vs_truth' + suffix + '.pdf'), dpi=190)
                if not full:
                    pdf.savefig(fig, dpi=190)
                plt.close(fig)
    (OUT / 'beta_truth_scatter_audit.json').write_text(json.dumps(audit, indent=2)+'\n')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
