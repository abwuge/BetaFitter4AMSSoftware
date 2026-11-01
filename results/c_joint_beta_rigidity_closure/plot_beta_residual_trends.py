#!/usr/bin/env python3
"""Fit beta residual cores in true-beta bins and plot mean and sigma."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy.optimize import minimize
from scipy.special import ndtr

HERE = Path(__file__).resolve().parent / 'expanded_validation'
OUT = HERE / 'figures'
EDGES = np.linspace(.68, .98, 16)
METHODS = [('beta_tof', 'u_raw', r'$\beta_{\rm TOF}$ (linear)', '#5790fc'),
           ('beta_masspeak', 'u_joint_anchored', r'$\beta_{\rm masspeak}$ (anchored joint)', '#f89c20')]


def fit_core(residual):
    """Conditional Gaussian MLE in a fixed median +/- 2.5 W68 window."""
    lo, median, hi = np.percentile(residual, [16, 50, 84])
    width = (hi-lo)/2
    left, right = median-2.5*width, median+2.5*width
    x = residual[(residual >= left) & (residual <= right)]

    def nll(p):
        mu, logsigma = p
        sigma = np.exp(logsigma)
        norm = ndtr((right-mu)/sigma)-ndtr((left-mu)/sigma)
        if norm <= 0:
            return float('inf')
        return len(x)*(logsigma+np.log(norm)) + .5*np.sum(((x-mu)/sigma)**2)

    result = minimize(nll, [median, np.log(width)], method='L-BFGS-B',
                      bounds=[(left, right), (np.log(width/5), np.log(width*5))],
                      options={'ftol': 1e-12, 'gtol': 1e-7, 'maxiter': 1000})
    if not result.success:
        raise RuntimeError(result.message)
    p = result.x
    steps = np.array([max(np.exp(p[1])*1e-3, 1e-4), 1e-3])
    h = np.zeros((2, 2))
    for i in range(2):
        ei = np.eye(2)[i]*steps[i]
        h[i,i] = (nll(p+ei)-2*nll(p)+nll(p-ei))/steps[i]**2
        for j in range(i):
            ej = np.eye(2)[j]*steps[j]
            h[i,j] = h[j,i] = (nll(p+ei+ej)-nll(p+ei-ej)-nll(p-ei+ej)+nll(p-ei-ej))/(4*steps[i]*steps[j])
    assert np.all(np.linalg.eigvalsh(h)>0), 'Non-positive fit information'
    errors = np.sqrt(np.diag(np.linalg.inv(h)))
    return dict(mean=p[0], mean_error=errors[0], sigma=np.exp(p[1]),
                sigma_error=np.exp(p[1])*errors[1], n_core=len(x),
                fit_low=left, fit_high=right)


def main():
    with np.load(HERE/'predictions.npz') as f:
        data = dict(f)
    good = data['truth_valid'] & np.isfinite(data['u_true']) & (data['u_true']>0)
    truth = 1/data['u_true'][good]
    rows = []
    for name, key, label, color in METHODS:
        residual = 1000*(1/data[key][good]-truth)
        assert np.all(np.isfinite(residual))
        for low, high in zip(EDGES[:-1], EDGES[1:]):
            values = residual[(truth>=low) & (truth<high)]
            assert len(values)>=100
            fit = fit_core(values)
            rows.append(dict(method=name, beta_low=low, beta_high=high,
                             beta_center=(low+high)/2, n=len(values),
                             **fit, sample_mean=float(np.mean(values)),
                             sample_std=float(np.std(values,ddof=1))))
    with (HERE/'beta_residual_trends.tsv').open('w') as f:
        writer=csv.DictWriter(f, fieldnames=list(rows[0]),delimiter='\t')
        writer.writeheader();writer.writerows(rows)
    plt.rcParams.update({'font.size':13, 'axes.labelsize':16, 'axes.titlesize':16,
                         'legend.fontsize':12, 'xtick.direction':'in', 'ytick.direction':'in',
                         'xtick.top':True, 'ytick.right':True, 'pdf.fonttype':42})
    with PdfPages(OUT/'beta_residual_mean_sigma.pdf') as pdf:
        for metric, ylabel in [('mean',r'$\mu(\beta-\beta_{\rm truth})\;[10^{-3}]$'),
                               ('sigma',r'$\sigma(\beta-\beta_{\rm truth})\;[10^{-3}]$')]:
            fig, ax=plt.subplots(figsize=(9,6.5))
            for name,key,label,color in METHODS:
                selected=[r for r in rows if r['method']==name]
                ax.errorbar([r['beta_center'] for r in selected], [r[metric] for r in selected],
                            yerr=[r[metric+'_error'] for r in selected],
                            fmt='o-', markersize=4.5, lw=1.5, capsize=3, color=color, label=label)
            if metric=='mean':
                ax.axhline(0,color='0.5',lw=1,ls='--',zorder=0)
            ax.set(xlabel=r'$\beta_{\rm truth}$',ylabel=ylabel,xlim=(.68,.98),
                   title='C12 + C13: beta residual '+metric)
            ax.margins(y=.3)
            if metric=='sigma':
                ax.set_ylim(0,ax.get_ylim()[1])
            legend=ax.legend(frameon=False,loc='upper left')
            for text, (_,_,_,color) in zip(legend.get_texts(),METHODS): text.set_color(color)
            fig.tight_layout();pdf.savefig(fig)
            fig.savefig(OUT/('beta_residual_'+metric+'.png'),dpi=180)
            fig.savefig(OUT/('beta_residual_'+metric+'.pdf'))
            plt.close(fig)
    (HERE/'beta_residual_trends.json').write_text(json.dumps({
        'residual':'beta_reco - beta_truth', 'table_units':'1e-3 for residual values and fit errors',
        'beta_masspeak':'1 / u_joint_anchored', 'beta_tof':'1 / u_raw',
        'fit':'Unbinned conditional Gaussian likelihood in median +/- 2.5 W68; truncation normalization included',
        'errors':'1-sigma inverse observed-Hessian fit errors; no source-file resampling',
        'beta_edges':EDGES.tolist(), 'input_events':int(good.sum()),
        'events_in_plotted_bins':int(sum(r['n'] for r in rows if r['method']=='beta_tof')),
        'sample_moments':'Untrimmed sample_mean and sample_std also saved in TSV'
    },indent=2)+'\n')
    print('Saved 2 plots and %d fitted bin/method rows'%len(rows))


if __name__=='__main__':
    main()
