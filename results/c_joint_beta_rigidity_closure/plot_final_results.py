#!/usr/bin/env python3
"""Plot saved, frozen C12/C13 validation predictions without refitting."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'figures'
COLORS = ['#5790fc', '#f89c20', '#e42536', '#964a8b']
LABELS = ['Raw', r'Frozen $\beta$-only', 'Unanchored joint', 'Anchored joint']
MASSES = np.array([11.17793, 12.11221])


def main():
    OUT.mkdir(exist_ok=True)
    with np.load(HERE.parent / 'c_isotope_proportion_closure/controlled_cache.npz') as f:
        d = {k[7:]: f[k] for k in f.files if k.startswith('valid__')}
    with np.load(HERE / 'joint_identifiability_predictions.npz') as f:
        a = dict(f)
    with np.load(HERE / 'anchored_joint_predictions.npz') as f:
        b = dict(f)
    for source in (a, b):
        for key in ('run', 'event', 'label', 'truth_valid'):
            assert np.array_equal(source[key], d[key]), key
    valid = d['truth_valid'] & np.isfinite(d['u_true']) & np.isfinite(d['R_true']) & (d['R_true'] > 0)
    assert valid.sum() == 6013
    pairs = [(a['u_raw'], a['R_raw']), (a['u_beta_only_frozen'], a['R_beta_only_frozen']),
             (a['u_joint_balanced'], a['R_joint_balanced']),
             (b['u_anchored_joint_balanced'], b['R_anchored_joint_balanced'])]
    values = []
    rows = []
    for name, (u, r) in zip(LABELS, pairs):
        m2 = (6*r)**2 * (u*u-1)
        m = np.full(len(u), np.nan)
        m[m2 >= 0] = np.sqrt(m2[m2 >= 0])
        obs = {'mass': m, 'mass2': m2, 'du': u-d['u_true'],
               'dbeta': 1/u-1/d['u_true'], 'dlogR': np.log(r/d['R_true'])}
        values.append(obs)
        for sample, mask in [('all', valid), ('C12', valid & (d['label']==0)), ('C13', valid & (d['label']==1))]:
            for key, arr in obs.items():
                x = arr[mask & np.isfinite(arr)]
                q = np.percentile(x, [16, 50, 84])
                rows.append(dict(method=name, sample=sample, observable=key, n=len(x), median=q[1], w68=(q[2]-q[0])/2))
    with (OUT/'summary.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    # Check the saved evaluation convention against the original numerical report.
    ref = json.loads((HERE/'joint_cluster_bootstrap_results.json').read_text())['observed']
    for name, obs in zip(['raw','beta_only_frozen','joint_unanchored','joint_anchored'], values):
        for key, metric in [('du','w68_du'), ('dlogR','w68_dlogR')]:
            q = np.percentile(obs[key][valid], [16,84])
            assert np.isclose((q[1]-q[0])/2, ref[name][metric], rtol=1e-10)
    plt.rcParams.update({'font.size': 13, 'axes.labelsize': 15, 'axes.titlesize': 16,
                         'legend.fontsize': 11, 'xtick.direction':'in', 'ytick.direction':'in',
                         'xtick.top':True, 'ytick.right':True, 'pdf.fonttype':42})
    audit = []
    with PdfPages(OUT/'final_results.pdf') as pdf:
        def draw(filename, key, mask, bins, xlabel, title, scale=1):
            fig, ax = plt.subplots(figsize=(9,6.5))
            for i, obs in enumerate(values):
                x = obs[key][mask]*scale
                x = x[np.isfinite(x)]
                q = np.percentile(x,[16,50,84]); width=(q[2]-q[0])/2
                ax.hist(x, bins=bins, histtype='step', lw=1.9, color=COLORS[i],
                        linestyle=['-','--','-.','-'][i],
                        label=f'{LABELS[i]}: median {q[1]:.4g}, $W_{{68}}$ {width:.4g}')
                audit.append(dict(figure=filename,method=LABELS[i],n=len(x),outside=int(np.sum((x<bins[0])|(x>bins[-1])))))
            ax.set(xlabel=xlabel, ylabel=f'Events / {bins[1]-bins[0]:.3g}', title=title, xlim=(bins[0],bins[-1]))
            ax.set_ylim(0, ax.get_ylim()[1]*1.36)
            legend = ax.legend(loc='upper right',frameon=False)
            for text, color in zip(legend.get_texts(), COLORS):
                text.set_color(color)
            ax.text(.025,.96,f'N = {mask.sum():,}',transform=ax.transAxes,va='top',fontsize=11)
            fig.tight_layout(); pdf.savefig(fig); fig.savefig(OUT/(filename+'.png'),dpi=170); plt.close(fig)
        for isotope in (0,1):
            mask=valid & (d['label']==isotope)
            draw(f'mass_C{12+isotope}','mass',mask,np.arange(4,20.0001,.25),
                 r'Reconstructed mass [GeV/$c^2$]',rf'$^{{{12+isotope}}}$C mass spectrum',)
        # Equal-area shapes compare isotope separation independently of sample sizes.
        fig, ax=plt.subplots(figsize=(9,6.5))
        bins=np.arange(4,20.0001,.25)
        for i,(method,isotope) in enumerate([(0,0),(0,1),(3,0),(3,1)]):
            mask=valid & (d['label']==isotope)
            x=values[method]['mass'][mask]; x=x[np.isfinite(x)]
            ax.hist(x,bins=bins,weights=np.ones(len(x))/len(x)/(bins[1]-bins[0]),
                    histtype='step',lw=2,color=COLORS[i],linestyle='--' if method==0 else '-',
                    label=rf'$^{{{12+isotope}}}$C, '+('Raw' if method==0 else 'Anchored joint'))
        ax.set(xlabel=r'Reconstructed mass [GeV/$c^2$]',ylabel=r'Probability density [(GeV/$c^2$)$^{-1}$]',
               title='C12 / C13 mass separation',xlim=(4,20))
        ax.set_ylim(0,ax.get_ylim()[1]*1.24)
        legend = ax.legend(frameon=False,loc='upper right')
        for text, color in zip(legend.get_texts(), COLORS):
            text.set_color(color)
        fig.tight_layout(); pdf.savefig(fig); fig.savefig(OUT/'mass_separation.png',dpi=170); plt.close(fig)
        for sample,mask in [('all',valid),('C13',valid & (d['label']==1))]:
            draw('inverse_beta_'+sample,'du',mask,np.linspace(-.05,.05,81)*1000,
                 r'$10^3(1/\beta_{\rm reco}-1/\beta_{\rm true})$',
                 'Inverse-beta residual: '+('C12 + C13' if sample=='all' else 'C13'),1000)
        draw('beta_C13','dbeta',valid & (d['label']==1),np.linspace(-.03,.03,81)*1000,
             r'$10^3(\beta_{\rm reco}-\beta_{\rm true})$','Beta residual: C13',1000)
        draw('rigidity_all','dlogR',valid,np.linspace(-.5,.5,81),
             r'$\ln(R_{\rm reco}/R_{\rm true})$','Rigidity residual: C12 + C13')
    (OUT/'plot_audit.json').write_text(json.dumps({'events':int(valid.sum()),
        'C12':int(np.sum(valid & (d['label']==0))), 'C13':int(np.sum(valid & (d['label']==1))),
        'w68_definition':'(q84-q16)/2; computed over full sample, including tails outside displayed range',
        'histograms':audit},indent=2)+'\n')
    print(OUT/'final_results.pdf')


if __name__ == '__main__':
    main()
