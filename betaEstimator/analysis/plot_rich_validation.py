#!/usr/bin/env python3
"""Numerical validation plots, one old/new comparison per figure."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import rich_transport as rt
import evaluate_rich_transport as ev

DEST=rt.ROOT/'results/rich_transport'
FIG=DEST/'figures'


def load(name):
    with np.load(str(DEST/name)) as f:return {k:f[k] for k in f.files}


def histogram(name,title,old,new,xlabel,scale=1000.,limit=40.,step=1.):
    old=old*scale;new=new*scale
    edges=np.arange(-limit,limit+step*.5,step)
    fig,ax=plt.subplots(figsize=(9.2,6.4))
    for value,label in [(old,'Linear TOF'),(new,'Anchored energy')]:
        text='%s\nW68 = %.3f; median = %+.3f'%(label,ev.width(value),np.median(value))
        ax.hist(value,bins=edges,histtype='step',linewidth=1.8,label=text)
    ax.set_xlim(-limit,limit);ax.set_ylim(bottom=0,top=ax.get_ylim()[1]*1.23)
    ax.set_xlabel(xlabel,fontsize=15);ax.set_ylabel('Events / %.3g'%(step/scale),fontsize=15)
    ax.set_title(title,fontsize=16,pad=13)
    ax.tick_params(labelsize=12)
    ax.legend(loc='upper center',ncol=2,frameon=False,fontsize=11)
    ax.text(.98,.68,'N = %s'%format(len(old),','),ha='right',va='top',transform=ax.transAxes,fontsize=12)
    fig.tight_layout();fig.savefig(str(FIG/(name+'.png')),dpi=170);plt.close(fig)


def trend(name,title,p,mask,reference):
    edges=np.array([.78,.82,.86,.90,.94,.97,1.02])
    fig,ax=plt.subplots(figsize=(9.2,6.4))
    for key,label in [('raw_beta_planes','Linear TOF + transport'),('new_beta_planes','Anchored energy + transport')]:
        xs,ys=[],[]
        for low,high in zip(edges[:-1],edges[1:]):
            s=mask&(reference>=low)&(reference<high)
            if s.sum()>=30:
                xs.append(np.median(reference[s]));ys.append(1000*ev.width(p[key][s,-1]-reference[s]))
        ax.plot(xs,ys,marker='o',linewidth=1.7,label=label)
    ax.set_xlabel(r'$\beta_{\rm RICH}$',fontsize=15);ax.set_ylabel(r'$W_{68}(\hat\beta_{\rm RICH}-\beta_{\rm RICH})\times10^3$',fontsize=15)
    ax.set_title(title,fontsize=16,pad=13);ax.tick_params(labelsize=12);ax.legend(frameon=False,fontsize=12)
    fig.tight_layout();fig.savefig(str(FIG/(name+'.png')),dpi=170);plt.close(fig)


def main():
    FIG.mkdir(exist_ok=True)
    mc=load('mc_predictions.npz');mq=load('mc_reference_v2_masks.npz')
    data=load('data_predictions.npz');dq=load('data_reference_v2_masks.npz')
    m=mc['truth_valid']&np.isfinite(mc['u_true'])&(mc['u_true']>0)
    histogram('mc_center_all','MC center reference (L5)',1/mc['u_raw'][m]-1/mc['u_true'][m],1/mc['u_new'][m]-1/mc['u_true'][m],r'$(\hat\beta_{L5}-\beta_{L5}^{\rm true})\times10^3$')
    matched=m&mq['charge_consistent_reference']
    histogram('mc_center_matched','MC center truth: charge-consistent RICH events',1/mc['u_raw'][matched]-1/mc['u_true'][matched],1/mc['u_new'][matched]-1/mc['u_true'][matched],r'$(\hat\beta_{L5}-\beta_{L5}^{\rm true})\times10^3$')
    for kind,p,q in [('mc',mc,mq),('data',data,dq)]:
        sel=q['charge_consistent_reference']
        histogram(kind+'_rich_charge_consistent',('MC' if kind=='mc' else 'Flight holdout')+': RICH reference, 4 < Q < 8',p['raw_beta_planes'][sel,-1]-p['rich_beta'][sel],p['new_beta_planes'][sel,-1]-p['rich_beta'][sel],r'$(\hat\beta_{\rm RICH}-\beta_{\rm RICH})\times10^3$',step=1. if kind=='mc' else 2.)
        trend(kind+'_rich_width_by_beta',('MC' if kind=='mc' else 'Flight holdout')+': charge-consistent RICH',p,sel,p['rich_beta'])
    # The reference audit shows all recorded events; the displayed interval
    # only sets the plot viewport, and overflow fractions remain in JSON.
    s=np.isfinite(mq['exact_center_rich_beta'])
    fig,ax=plt.subplots(figsize=(9.2,6.4))
    r=(mq['exact_center_rich_beta'][s]-mc['rich_beta'][s])*1000
    ax.hist(r,bins=np.linspace(-5,5,101),histtype='step',linewidth=1.8)
    ax.set_title('MC transport from exact center: RICH reference check',fontsize=16,pad=13)
    ax.set_xlabel(r'$(\beta_{\rm RICH}^{\rm exact\ center\ transport}-\beta_{\rm RICH})\times10^3$',fontsize=14)
    ax.set_ylabel('Events / 0.0001',fontsize=15);ax.tick_params(labelsize=12)
    ax.text(.97,.94,'N = %s\nW68 = %.3f\nmedian = %+.3f\nOutside window: %.2f%%'%(format(len(r),','),ev.width(r),np.median(r),100*np.mean(np.abs(r)>5)),ha='right',va='top',transform=ax.transAxes,fontsize=12)
    fig.tight_layout();fig.savefig(str(FIG/'mc_reference_core_check.png'),dpi=170);plt.close(fig)
    # Each point compares widths on exactly the same events before/after.
    allmc=json.loads((DEST/'mc_results.json').read_text())
    mv=json.loads((DEST/'mc_reference_v2_results.json').read_text())
    dv=json.loads((DEST/'data_reference_v2_results.json').read_text())
    rows=[('MC: all center truth',allmc['mc_center_all'],allmc['mc_center_bootstrap']),
          ('MC: matched center truth',mv['matched_center_truth'],mv['matched_center_truth']['bootstrap']),
          ('MC: matched RICH',mv['matched_RICH'],mv['matched_RICH']['bootstrap']),
          ('Flight holdout: RICH',dv['rich_all'],dv['rich_all']['bootstrap'])]
    fig,ax=plt.subplots(figsize=(9.6,5.2))
    for i,(label,result,boot) in enumerate(rows):
        ci=boot['beta_gain_pct_q025_q50_q975'];point=result['gain_beta_pct']
        # An asymmetric interval can exclude the point in a tiny sample;
        # draw its endpoints directly instead of forcing a negative xerr.
        ax.plot([ci[0],point,ci[2]],[i,i,i],linewidth=1.8,marker='o',markevery=[1],markersize=7)
    ax.set_yticks(range(len(rows)));ax.set_yticklabels([x[0] for x in rows],fontsize=12)
    ax.invert_yaxis();ax.set_xlabel('Reduction of beta residual W68 (%)',fontsize=14)
    ax.set_title('Paired width changes: source-cluster 95% intervals',fontsize=15,pad=12)
    ax.tick_params(axis='x',labelsize=12);fig.tight_layout();fig.savefig(str(FIG/'resolution_gain_summary.png'),dpi=170);plt.close(fig)
    print('FIGURES',len(list(FIG.glob('*.png'))),flush=True)

if __name__=='__main__':main()
