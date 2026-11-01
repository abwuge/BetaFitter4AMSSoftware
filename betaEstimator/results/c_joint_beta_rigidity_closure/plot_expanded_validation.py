#!/usr/bin/env python3
"""Plot additional-file validation and paired source-file uncertainty."""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from plot_final_results import COLORS, LABELS
from run_expanded_validation import METHODS, rc

HERE=Path(__file__).resolve().parent/'expanded_validation'
FIG=HERE/'figures'


def bootstrap(d, obs, nrep=500):
    """Resample whole files within isotope, using the same draws for all models."""
    ids=np.unique(np.column_stack([d['label'],d['source']]),axis=0)
    lookup={tuple(pair):i for i,pair in enumerate(ids)}
    cluster=np.array([lookup[(a,b)] for a,b in zip(d['label'],d['source'])])
    groups=[np.flatnonzero(ids[:,0]==a) for a in (0,1)]
    truth=d['truth_valid'] & np.isfinite(d['u_true']) & (d['R_true']>0)
    prepared={}
    for b in [-1,0,1,2,3]:
        sel=truth & ((d['bin']==b) if b>=0 else True)
        for m in range(4):
            for key in ['du','mass','dlogR','auc']:
                mask=sel & ((d['label']==1) if key in ['du','mass'] else True)
                x=obs[m]['mass2' if key=='auc' else key]
                mask &= np.isfinite(x)
                rows=np.flatnonzero(mask); rows=rows[np.argsort(x[rows],kind='stable')]
                starts=np.r_[0,np.flatnonzero(x[rows][1:]!=x[rows][:-1])+1]
                prepared[b,m,key]=(x[rows],cluster[rows],d['label'][rows],starts)
    results={k:[] for k in prepared}
    rng=np.random.RandomState(20260930)
    for rep in range(nrep):
        mult=np.zeros(len(ids))
        for g in groups:
            mult+=np.bincount(rng.choice(g,size=len(g),replace=True),minlength=len(ids))
        for key,(x,c,label,starts) in prepared.items():
            w=mult[c]
            if key[2]=='auc':
                neg=np.add.reduceat(w*(label==0),starts)
                pos=np.add.reduceat(w*(label==1),starts)
                val=np.sum(pos*(np.cumsum(neg)-.5*neg))/(pos.sum()*neg.sum())
            else:
                cum=np.cumsum(w)
                q=x[np.searchsorted(cum,np.array([.16,.84])*cum[-1])]
                val=(q[1]-q[0])/2
            results[key].append(val)
        if (rep+1)%100==0:
            print('Cluster bootstrap %d/%d'%(rep+1,nrep),flush=True)
    out=[]
    for b in [-1,0,1,2,3]:
        for m in range(1,4):
            for key in ['du','mass','dlogR','auc']:
                base=np.array(results[b,0,key]); val=np.array(results[b,m,key])
                delta=val-base if key=='auc' else 100*(1-val/base)
                q=np.percentile(delta,[2.5,50,97.5])
                out.append(dict(bin=b,method=METHODS[m],metric=key,low=q[0],median=q[1],high=q[2]))
    (HERE/'cluster_bootstrap.json').write_text(json.dumps({'replicates':nrep,'unit':'source file; stratified by isotope; paired across methods',
        'width_samples':'C13 for du and mass; C12+C13 for dlogR','interval':'percentile 95%; frozen models',
        'differences_vs_raw':out},indent=2)+'\n')
    return out


def main():
    FIG.mkdir(exist_ok=True)
    with np.load(HERE/'predictions.npz') as f:
        d=dict(f)
    valid=d['truth_valid'] & np.isfinite(d['u_true']) & np.isfinite(d['R_true']) & (d['R_true']>0)
    obs=[]
    for m in METHODS:
        u,r=d['u_'+m],d['R_'+m]
        mass2=(6*r)**2*(u*u-1)
        mass=np.full(len(u),np.nan); mass[mass2>=0]=np.sqrt(mass2[mass2>=0])
        obs.append(dict(mass=mass,mass2=mass2,du=u-d['u_true'],dbeta=1/u-1/d['u_true'],dlogR=np.log(r/d['R_true'])))
    decomposition=[]
    for b in range(4):
        mask=valid & (d['label']==1) & (d['bin']==b)
        for name,u,r in [('measured',d['u_raw'],d['R_raw']),
                         ('truth_beta',d['u_true'],d['R_raw']),
                         ('truth_rigidity',d['u_raw'],d['R_true'])]:
            error=(6*r[mask])**2*(u[mask]**2-1)-rc.MASSES[1]**2
            decomposition.append(dict(bin=b,scenario=name,n=int(mask.sum()),w68_mass2_error=rc.w68(error)))
    rc.write_tsv(HERE/'resolution_decomposition.tsv',decomposition)
    intervals=bootstrap(d,obs)
    plt.rcParams.update({'font.size':13,'axes.labelsize':15,'axes.titlesize':16,'legend.fontsize':11,
                         'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,'pdf.fonttype':42})
    pages=[]
    audit=[]
    with PdfPages(FIG/'expanded_results.pdf') as pdf:
        def finish(fig,ax,name,colors=COLORS,headroom=1.35):
            if headroom:
                ax.set_ylim(0,ax.get_ylim()[1]*headroom)
            legend=ax.legend(frameon=False,loc='upper right')
            for text,color in zip(legend.get_texts(),colors): text.set_color(color)
            fig.tight_layout();pdf.savefig(fig);fig.savefig(FIG/(name+'.png'),dpi=170);plt.close(fig);pages.append(name)
        def hist(key,mask,bins,xlabel,title,name,scale=1):
            fig,ax=plt.subplots(figsize=(9,6.5))
            for i in range(4):
                x=obs[i][key][mask]*scale;x=x[np.isfinite(x)]
                q=np.percentile(x,[16,50,84])
                ax.hist(x,bins=bins,histtype='step',lw=1.9,color=COLORS[i],linestyle=['-','--','-.','-'][i],
                        label=f'{LABELS[i]}: median {q[1]:.4g}, $W_{{68}}$ {(q[2]-q[0])/2:.4g}')
                audit.append(dict(figure=name,method=METHODS[i],n=len(x),outside=int(np.sum((x<bins[0])|(x>bins[-1])))))
            ax.set(xlabel=xlabel,ylabel=f'Events / {bins[1]-bins[0]:.3g}',title=title,xlim=(bins[0],bins[-1]))
            ax.text(.025,.96,f'N = {mask.sum():,}',transform=ax.transAxes,va='top',fontsize=11)
            finish(fig,ax,name)
        for b in [-1,0,1,2,3]:
            sel=valid & ((d['bin']==b) if b>=0 else True)
            rlabel='2.5–6 GV' if b<0 else f'{rc.R_EDGES[b]:g}–{rc.R_EDGES[b+1]:g} GV'
            suffix='all' if b<0 else 'bin%d'%b
            fig,ax=plt.subplots(figsize=(9,6.5)); bins=np.arange(4,20.0001,.2)
            for i,(m,label) in enumerate([(0,0),(0,1),(3,0),(3,1)]):
                x=obs[m]['mass'][sel & (d['label']==label)]; x=x[np.isfinite(x)]
                ax.hist(x,bins=bins,weights=np.ones(len(x))/(len(x)*.2),histtype='step',lw=2,
                        color=COLORS[i],linestyle='--' if m==0 else '-',
                        label=rf'$^{{{12+label}}}$C, '+('Raw' if m==0 else 'Anchored joint'))
            ax.set(xlabel=r'Reconstructed mass [GeV/$c^2$]',ylabel=r'Probability density [(GeV/$c^2$)$^{-1}$]',
                   title='C12 / C13 mass separation: '+rlabel,xlim=(4,20))
            ax.text(.025,.96,f'C12: {np.sum(sel & (d["label"]==0)):,}\nC13: {np.sum(sel & (d["label"]==1)):,}',transform=ax.transAxes,va='top',fontsize=11)
            finish(fig,ax,'mass_separation_'+suffix,headroom=1.25)
            hist('du',sel & (d['label']==1),np.linspace(-50,50,101),r'$10^3(1/\beta_{\rm reco}-1/\beta_{\rm true})$',
                 'C13 inverse-beta residual: '+rlabel,'inverse_beta_C13_'+suffix,1000)
        hist('mass',valid & (d['label']==1),np.arange(4,20.0001,.2),r'Reconstructed mass [GeV/$c^2$]','C13 mass spectrum: 2.5–6 GV','mass_C13')
        hist('dbeta',valid & (d['label']==1),np.linspace(-30,30,101),r'$10^3(\beta_{\rm reco}-\beta_{\rm true})$','C13 beta residual: 2.5–6 GV','beta_C13',1000)
        hist('dlogR',valid,np.linspace(-.5,.5,101),r'$\ln(R_{\rm reco}/R_{\rm true})$','C12 + C13 rigidity residual: 2.5–6 GV','rigidity_all')
        for key,title,ylabel in [('du','C13 inverse-beta resolution',r'$W_{68}$ reduction vs. raw [%]'),
                                  ('mass','C13 mass resolution',r'$W_{68}$ reduction vs. raw [%]'),
                                  ('auc','C12 / C13 mass discrimination',r'Mass$^2$ AUC increase vs. raw')]:
            fig,ax=plt.subplots(figsize=(9,6.5))
            for m in [1,2,3]:
                rows=[next(r for r in intervals if r['metric']==key and r['bin']==b and r['method']==METHODS[m]) for b in range(4)]
                y=np.array([r['median'] for r in rows]);low=np.array([r['low'] for r in rows]);high=np.array([r['high'] for r in rows])
                ax.errorbar(np.arange(4)+(m-2)*.12,y,yerr=[y-low,high-y],fmt='o',capsize=4,color=COLORS[m],label=LABELS[m])
            ax.axhline(0,color='0.6',lw=1,zorder=0)
            ax.set_xticks(np.arange(4));ax.set_xticklabels(['2.5–3','3–4','4–5','5–6'])
            ax.set(xlabel='Reconstructed rigidity [GV]',ylabel=ylabel,title=title)
            ax.margins(y=.3)
            ax.text(.025,.025,'95% source-file bootstrap intervals',transform=ax.transAxes,fontsize=11)
            finish(fig,ax,key+'_improvement',colors=COLORS[1:],headroom=0)
    (FIG/'audit.json').write_text(json.dumps({'pages':pages,'histograms':audit},indent=2)+'\n')
    print('Completed %d pages: %s'%(len(pages),FIG/'expanded_results.pdf'),flush=True)


if __name__=='__main__':
    main()
