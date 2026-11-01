#!/usr/bin/env python3
"""Compare frozen tracker predictions on the untouched whole-file test split."""
import csv
import json
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'c_joint_beta_rigidity_closure'))
from plot_beta_residual_trends import fit_core
from summarize_strict import load_stage
OUT=HERE/'figures'
COLORS=['#5790fc','#f89c20','#e42536','#964a8b']
EDGES=np.linspace(.68,.98,16)
MASSES=np.array([11.17793,12.11221])

def main():
    with np.load(str(HERE.parent/'c_mass_toy_oracle'/'predictions.npz')) as z:
        idx=np.flatnonzero(z['role']==b'test')
        data={k:z[k][idx] for k in ['label','source','u','u_true','R','R_true']}
    predictions=load_stage('strict_geometry')
    mapping={int(i):float(v) for i,v in zip(predictions['event_index'],predictions['BDTG'])}
    correction=np.array([mapping[int(i)] for i in idx])
    assert len(mapping)==len(idx) and np.all(np.isfinite(correction))
    beta=1/data['u_true']
    rigidity=[data['R'],data['R']*np.exp(-correction)]
    residuals={
      'rigidity':[100*np.log(r/data['R_true']) for r in rigidity],
      'mass2':[(6*r)**2*(data['u']**2-1)-MASSES[data['label']]**2 for r in rigidity]}
    series=[(0,0,'C12 raw'),(0,1,'C12 corrected'),(1,0,'C13 raw'),(1,1,'C13 corrected')]
    rows=[]
    for observable,values in residuals.items():
        for lab,method,label in series:
            for low,high in zip(EDGES[:-1],EDGES[1:]):
                mask=(data['label']==lab)&(beta>=low)&(beta<high)
                if mask.sum()<80: continue
                fit=fit_core(values[method][mask])
                rows.append(dict(observable=observable,isotope=lab,method=method,
                  beta_low=low,beta_high=high,beta_center=(low+high)/2,n=int(mask.sum()),**fit))
    with (HERE/'residual_trends.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    plt.rcParams.update({'font.size':13,'axes.labelsize':16,'axes.titlesize':16,
       'legend.fontsize':12,'xtick.direction':'in','ytick.direction':'in',
       'xtick.top':True,'ytick.right':True,'pdf.fonttype':42})
    with PdfPages(str(OUT/'current_residual_results.pdf')) as pdf:
        for observable in ['rigidity','mass2']:
            for metric in ['mean','sigma']:
                fig,ax=plt.subplots(figsize=(9,6.5))
                for color,(lab,method,label) in zip(COLORS,series):
                    rr=[r for r in rows if r['observable']==observable and r['isotope']==lab and r['method']==method]
                    ax.errorbar([r['beta_center'] for r in rr],[r[metric] for r in rr],
                      yerr=[r[metric+'_error'] for r in rr],fmt=('o-' if lab==0 else 's--'),
                      color=color,lw=1.8,markersize=4.5,capsize=3,label=label)
                symbol=r'\mu' if metric=='mean' else r'\sigma'
                ylabel=(r'$'+symbol+r'[\ln(R/R_{\rm truth})]\;[\%]$' if observable=='rigidity'
                  else r'$'+symbol+r'(m^2-m^2_{\rm truth})\;[\mathrm{GeV}^2/c^4]$')
                ax.set(xlabel=r'$\beta_{\rm truth}$',ylabel=ylabel,xlim=(.68,.98),
                       title=('Tracker rigidity' if observable=='rigidity' else 'Mass squared, measured TOF beta')+' residual '+metric)
                ax.margins(y=.3)
                if metric=='sigma':ax.set_ylim(0,ax.get_ylim()[1])
                legend=ax.legend(frameon=False,ncol=2,loc='upper left')
                for text,color in zip(legend.get_texts(),COLORS):text.set_color(color)
                fig.tight_layout();pdf.savefig(fig)
                stem=observable+'_residual_'+metric
                fig.savefig(str(OUT/(stem+'.png')),dpi=180);fig.savefig(str(OUT/(stem+'.pdf')))
                plt.close(fig)
        fig,ax=plt.subplots(figsize=(9,6.5))
        bins=np.arange(40,260.01,2.)
        for color,(lab,method,label) in zip(COLORS,series):
            m2=(6*rigidity[method])**2*(data['u']**2-1)
            ax.hist(m2[data['label']==lab],bins=bins,histtype='step',color=color,
              lw=1.8,linestyle=('-' if lab==0 else '--'),label=label)
        ax.set(xlabel=r'$m^2\;[\mathrm{GeV}^2/c^4]$',ylabel=r'Events / $(2\;\mathrm{GeV}^2/c^4)$',
          title='C12 + C13 mass peaks: raw and corrected tracker rigidity',xlim=(40,260))
        legend=ax.legend(frameon=False,ncol=2)
        for text,color in zip(legend.get_texts(),COLORS):text.set_color(color)
        fig.tight_layout();pdf.savefig(fig);fig.savefig(str(OUT/'mass2_peaks.png'),dpi=180)
        plt.close(fig)
    (HERE/'residual_trends.json').write_text(json.dumps({
      'test_events':len(idx),'C12':int(sum(data['label']==0)),'C13':int(sum(data['label']==1)),
      'model':'strict_geometry; MC-truth-supervised tracker oracle',
      'beta_edges':EDGES.tolist(),'truth_reference':'Tracker L5, mevmom1[10]',
      'fit':'Conditional Gaussian MLE, median +/- 2.5 W68 window',
      'errors':'Observed-Hessian 1-sigma fit errors',
      'rigidity_residual':'100*ln(R/R_true), percent',
      'mass2_residual':'(6*R)^2*(u_raw^2-1)-isotope_mass^2, GeV^2/c^4',
      'mass_histogram':'raw counts, same test events for each method','pages':5},indent=2)+'\n')
    print('Saved 5 pages and',len(rows),'fit rows')
if __name__=='__main__': main()
