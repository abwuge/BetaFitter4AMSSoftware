#!/usr/bin/env python3
"""Plot the paired residual distributions underlying every wide-range point."""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import flight_beta_wide as w

COLORS = ['#5790fc', '#f89c20']
RESULT_DEST = w.DEST
OUT = RESULT_DEST / 'distributions'
COLLECTIONS_ONLY = False


def number(x):
    return ('%g' % x).replace('.', 'p')


def draw(axis, row, residuals, edges, full):
    fig, ax = plt.subplots(figsize=(8.8, 6.2))
    counts = []
    for values, name, color in zip(residuals, ['Linear TOF', 'Flight mass fit'], COLORS):
        bins, _ = np.histogram(values, bins=edges)
        outside = int(np.sum((values < edges[0]) | (values > edges[-1])))
        if len(values) >= 2:
            label = '%s\n$W_{68}$ = %.3f; mean = %+.3f' % (name, w.ev.width(values), np.mean(values))
        else:
            label = '%s\nInsufficient statistics (N = %d)' % (name, len(values))
        if not full:
            label += '\nOutside view: %d' % outside
        ax.hist(values, bins=edges, histtype='step', color=color, linewidth=1.8, label=label)
        counts.append(dict(histogram=bins, outside_view=outside))
        if full and outside:
            raise RuntimeError('Full-range histogram omitted events')
    ax.axvline(0, color='0.6', linestyle='--', linewidth=1., ymax=.8)
    ax.set_xlim(edges[0], edges[-1])
    peak = max(1, max(np.max(c['histogram']) for c in counts))
    if full:
        ax.set_ylim(.7, peak * 25)
        ax.set_yscale('log')
    else:
        ax.set_ylim(0, peak * 1.55)
    ax.set_xlabel(r'$(\widehat\beta_{\rm RICH}-\beta_{\rm RICH}^{\rm meas})\times10^3$', fontsize=15)
    ax.set_ylabel('Events / %.3f' % ((edges[1] - edges[0]) / 1000), fontsize=14)
    interval = (r'$%.3g\leq R<%.3g$ GV' if axis == 'rigidity' else
                r'$%.5g\leq\beta_{\rm RICH}<%.5g$') % (row['low'], row['high'])
    ax.set_title(interval + '; N = %d' % row['n'], fontsize=14, pad=12)
    ax.tick_params(labelsize=12)
    ax.legend(frameon=False, fontsize=10.5, loc='upper center', ncol=2)
    fig.tight_layout()
    return fig, counts


def main():
    w.check_freeze()
    OUT.mkdir(exist_ok=True)
    path = RESULT_DEST / 'predictions.npz'
    before = w.ev.digest(path)
    result_path = RESULT_DEST / 'axis_results.json'
    if not result_path.exists():
        result_path = RESULT_DEST / 'results.json'
    result = json.loads(result_path.read_text())
    with np.load(str(path), allow_pickle=False) as f:
        p = {k: f[k] for k in f.files}
    raw = (p['raw_beta_rich'] - p['rich_beta']) * 1000
    new = (p['new_beta_rich'] - p['rich_beta']) * 1000
    mask = p['reference_mask']
    all_values = np.r_[raw[mask], new[mask]]
    if not len(all_values):
        raise RuntimeError('No reference events to plot')
    core_edges = np.arange(-60, 60.1, 2.)
    global_extent = [float(all_values.min()), float(all_values.max())]
    records = []
    for axis, rows in result['rows'].items():
        values = p['R'] if axis == 'rigidity' else p['rich_beta']
        with PdfPages(str(OUT / (axis + '_residuals.pdf'))) as core_pdf, \
                PdfPages(str(OUT / (axis + '_residuals_full.pdf'))) as full_pdf:
            for row in rows:
                select = mask & (values >= row['low']) & (values < row['high'])
                if int(select.sum()) != row['n']:
                    raise RuntimeError('Point event count changed')
                residuals = [raw[select], new[select]]
                bin_values = np.concatenate(residuals)
                if len(bin_values):
                    lower = np.floor(bin_values.min() / 5) * 5
                    upper = max(lower + 5, np.ceil(bin_values.max() / 5) * 5)
                    full_edges = np.arange(lower, upper + 5.1, 5.)
                else:
                    full_edges = np.array([-5., 0., 5.])
                for vals, method in zip(residuals, ['baseline_beta', 'new_beta']):
                    if row['n'] >= 2 and not np.isclose(w.ev.width(vals) / 1000, row[method]['w68'], atol=1e-15, rtol=0):
                        raise RuntimeError('Histogram population differs from plotted point')
                stem = '%s_%s_%s' % (axis, number(row['low']), number(row['high']))
                record = dict(axis=axis, low=row['low'], high=row['high'], n=row['n'], stem=stem,
                              full_edges=full_edges)
                for full, edges, pages in [(False, core_edges, core_pdf), (True, full_edges, full_pdf)]:
                    name = stem + ('_full' if full else '')
                    fig, counts = draw(axis, row, residuals, edges, full)
                    if not COLLECTIONS_ONLY:
                        fig.savefig(str(OUT / (name + '.pdf')))
                        fig.savefig(str(OUT / (name + '.png')), dpi=170)
                    pages.savefig(fig)
                    plt.close(fig)
                    record['full' if full else 'core'] = counts
                records.append(record)
        if sum(r['n'] for r in rows) != int(mask.sum()):
            raise RuntimeError('Incomplete axis coverage')
    if w.ev.digest(path) != before:
        raise RuntimeError('Predictions changed during plotting')
    w.ev.dump(OUT / 'histograms.json', dict(prediction_sha256=before,
              script_sha256=w.ev.digest(Path(__file__)), core_edges=core_edges,
              full_edges_per_bin=True, global_residual_extent=global_extent, residual_units='beta residual times 1000', records=records))
    print('Complete:', len(records), 'bins;', 2 * len(records), 'paired distributions', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', type=Path, default=w.DEST)
    parser.add_argument('--collections-only', action='store_true')
    args = parser.parse_args()
    COLLECTIONS_ONLY = args.collections_only
    RESULT_DEST = args.results
    OUT = RESULT_DEST / 'distributions'
    main()
