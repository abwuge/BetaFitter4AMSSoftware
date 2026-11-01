#!/usr/bin/env python3
"""Separate reconstructed-rigidity and independent-reference-beta figures."""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import flight_beta_wide as experiment
from evaluate_flight_beta_energy import paired_errors

COLORS = ['#5790fc', '#f89c20']
PLOT_DEST = experiment.DEST
CLUSTER_KEY = 'source'


def axes(axis):
    fig, ax = plt.subplots(figsize=(8.8, 6.2))
    if axis == 'rigidity':
        ax.set_xscale('log')
        ax.set_xlim(2.5, 200)
        ax.set_xticks([2.5, 5, 10, 20, 50, 100, 200])
        ax.set_xticklabels(['2.5', '5', '10', '20', '50', '100', '200'])
        ax.set_xlabel(r'Reconstructed rigidity $R$ (GV)', fontsize=15)
    else:
        ax.set_xlim(.78, 1.02)
        ax.set_xticks([.8, .85, .9, .95, .98, 1.])
        ax.set_xlabel(r'Reference $\beta_{\rm RICH}$', fontsize=15)
    ax.tick_params(labelsize=12)
    return fig, ax


def save(fig, name):
    folder = PLOT_DEST / 'figures'
    folder.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(str(folder / (name + '.pdf')))
    fig.savefig(str(folder / (name + '.png')), dpi=170)
    if name.startswith('beta_'):
        ax = fig.axes[0]
        ax.set_xlim(.975, 1.001)
        ax.set_xticks([.98, .985, .99, .995, 1.])
        ax.set_xticklabels(['0.980', '0.985', '0.990', '0.995', '1.000'])
        fig.savefig(str(folder / (name + '_zoom.pdf')))
        fig.savefig(str(folder / (name + '_zoom.png')), dpi=170)
    plt.close(fig)


def main():
    experiment.check_freeze()
    result = json.loads((PLOT_DEST/'results.json').read_text())
    with np.load(str(PLOT_DEST/'predictions.npz'),allow_pickle=False) as f:
        p = {k:f[k] for k in f.files}
    for axis, rows in result['rows'].items():
        values = p['R'] if axis=='rigidity' else p['rich_beta']
        usable = []
        for row in rows:
            if row['n']<20:
                continue
            sel = p['reference_mask'] & (values>=row['low']) & (values<row['high'])
            row['errors'] = paired_errors(p['raw_beta_rich'][sel],p['new_beta_rich'][sel],p['rich_beta'][sel],p[CLUSTER_KEY][sel])
            expected = np.array(row['bootstrap']['beta_gain_pct_q025_q50_q975'])[[0,2]]
            if not np.allclose(expected,row['errors']['gain_95ci'],atol=1e-9,rtol=0):
                raise RuntimeError('Bootstrap reproduction failed')
            usable.append(row)
        x = [r['x'] for r in usable]
        for metric,error,label in [
                ('w68','W68',r'$W_{68}(\widehat\beta_{\rm RICH}-\beta_{\rm RICH})\times10^3$'),
                ('mean','mean',r'$\langle\widehat\beta_{\rm RICH}-\beta_{\rm RICH}\rangle\times10^3$')]:
            fig,ax = axes(axis)
            for method,name,color,marker,prefix in [
                    ('baseline_beta','Linear TOF',COLORS[0],'o','old'),
                    ('new_beta','Flight mass fit',COLORS[1],'s','new')]:
                ax.errorbar(x,[1000*r[method][metric] for r in usable],
                            yerr=[1000*r['errors'][prefix+'_'+error+'_sd'] for r in usable],
                            marker=marker,color=color,capsize=3,linewidth=1.6,label=name)
            ax.set_ylabel(label,fontsize=14)
            ax.legend(frameon=False,fontsize=12,loc='upper left',
                      title='Holdout N = %d; bootstrap SD'%result['reference_n'],title_fontsize=11)
            if metric=='w68':
                ax.set_ylim(0,max(1000*r[m]['w68'] for r in usable for m in ('baseline_beta','new_beta'))*1.35)
            else:
                ax.axhline(0,color='0.6',linestyle='--',linewidth=1)
                ymin,ymax=ax.get_ylim()
                ax.set_ylim(ymin,ymax+.22*(ymax-ymin))
            save(fig,axis+('_width' if metric=='w68' else '_bias'))
        fig,ax = axes(axis)
        for r in usable:
            ax.plot([r['x']]*2,r['errors']['gain_95ci'],color=COLORS[1],linewidth=1.8)
            ax.plot(r['x'],r['gain_beta_pct'],'o',color=COLORS[1])
        ax.axhline(0,color=COLORS[0],linestyle='--',linewidth=1.2)
        ax.set_ylabel(r'Reduction of $\beta$ residual $W_{68}$ (%)',fontsize=14)
        ax.text(.03,.97,'N = %d\nPaired bootstrap: 95%% intervals'%result['reference_n'],
                transform=ax.transAxes,va='top',fontsize=11)
        ymin,ymax=ax.get_ylim()
        ax.set_ylim(ymin,ymax+.22*(ymax-ymin))
        save(fig,axis+'_gain')
    experiment.ev.dump(PLOT_DEST/'axis_results.json',result)
    print('SIX FIGURES + THREE BETA ZOOMS COMPLETE',flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', type=Path, default=experiment.DEST)
    parser.add_argument('--cluster', choices=['source', 'run'], default='source')
    args = parser.parse_args()
    PLOT_DEST = args.results
    CLUSTER_KEY = args.cluster
    main()
