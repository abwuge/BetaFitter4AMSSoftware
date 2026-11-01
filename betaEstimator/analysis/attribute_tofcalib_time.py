#!/usr/bin/env python3
"""Measured-covariate standardization of a frozen model's calendar performance.

All balancing uses measured covariates, not the beta residual or model gain.
Later timing/energy controls are mediation sensitivities, not causal proof.
"""
import argparse
from datetime import datetime, timezone
import json
import multiprocessing as mp
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import flight_beta_wide as w
from prepare_tofcalib_validation import OUT

P=C=YEAR=REF=DIMS=None
STAGES=['unadjusted','kinematics','reference_quality','geometry','energy_inputs','timing_inputs']
REPEATS=200


def quantile(values, weights, probabilities):
    positive=weights>0
    x=values[positive];v=weights[positive]
    if not len(x):raise RuntimeError('No weighted observations')
    order=np.argsort(x,kind='stable');x=x[order];v=v[order]
    if np.ptp(v)<1e-12:return np.percentile(x,np.array(probabilities)*100)
    return np.interp(probabilities,(np.cumsum(v)-.5*v)/v.sum(),x)


def measure(indices, weights):
    old=P['raw_beta_rich'][indices]-P['rich_beta'][indices]
    new=P['new_beta_rich'][indices]-P['rich_beta'][indices]
    widths=[float(np.diff(quantile(x,weights,[.16,.84]))[0]/2) for x in (old,new)]
    return dict(n=len(indices),effective_n=float(weights.sum()**2/np.sum(weights**2)),
        old_w68=widths[0],new_w68=widths[1],gain_pct=100*(1-widths[1]/widths[0]),
        old_bias=float(np.average(old,weights=weights)),new_bias=float(np.average(new,weights=weights)))


def category(value, ref):
    edges=np.unique(np.percentile(value[ref],[20,40,60,80]))
    return np.searchsorted(edges,value,side='right').astype(np.int16)


def build_dimensions():
    r=np.searchsorted(w.R_EDGES,P['R'],side='right')-1
    b=np.clip(np.searchsorted(w.BETA_EDGES,P['rich_beta'],side='right')-1,0,len(w.BETA_EDGES)-2)
    joint=(r*2+P['NaF'].astype(int))*(len(w.BETA_EDGES)-1)+b
    logedges=np.geomspace(2.5,200,25)
    k=[('R_beta_radiator_joint',joint.astype(np.int16)),
       ('rigidity_finer',np.searchsorted(logedges,P['R'],side='right').astype(np.int16)),
       ('beta_finer',np.searchsorted(np.array([.78,.82,.85,.9,.95,.955,.96,.965,.97,.975,.98,.985,.99,.9925,.995,.9965,.998,.999,1,1.001,1.002,1.02]),P['rich_beta'],side='right').astype(np.int16))]
    quality=[(key,category(C[key],REF)) for key in ['rich_pb','rich_error','rich_npe_expected','rich_npe_ratio',
       'rich_beta_consistency','rich_charge_consistency','rich_pmt','rich_used','rich_charge','relR','chi_x','chi_y']]
    cont=C['continuous']
    geometry=[('bar%d'%i,C['bars'][:,i].astype(np.int16)) for i in range(4)]
    geometry += [('geometry%d'%i,category(cont[:,i],REF)) for i in [0,3,4,7,8,9,29]]
    geometry += [(key,category(C[key],REF)) for key in ['rich_x','rich_y']]
    energy=[('energy%d'%i,category(cont[:,i],REF)) for i in range(18,29)]
    timing=[('timing%d'%i,category(cont[:,i],REF)) for i in range(14,18)]
    return [[],k,k+quality,k+quality+geometry,k+quality+geometry+energy,k+quality+geometry+energy+timing]


def supported(ref, target, dims):
    # Drop categories absent in either cohort; never select on outcomes.
    for _ in range(20):
        old=(len(ref),len(target))
        for _,values in dims:
            common=np.intersect1d(np.unique(values[ref]),np.unique(values[target]))
            ref=ref[np.isin(values[ref],common)];target=target[np.isin(values[target],common)]
        if old==(len(ref),len(target)):return ref,target
    raise RuntimeError('Common support failed to stabilize')


def balance(ref, target, dims):
    weights=np.ones(len(target),float)
    margins=[]
    for name,values in dims:
        n=int(values.max())+1
        a=np.bincount(values[ref],minlength=n)/len(ref)
        margins.append((name,values[target],a))
    error=0
    for iteration in range(300):
        for _,bins,desired in margins:
            observed=np.bincount(bins,weights=weights,minlength=len(desired))/weights.sum()
            ratio=np.divide(desired,observed,out=np.zeros_like(desired),where=observed>0)
            weights *= ratio[bins]
            weights *= len(weights)/weights.sum()
        error=max([np.max(np.abs(np.bincount(bins,weights=weights,minlength=len(desired))/weights.sum()-desired))
                   for _,bins,desired in margins]+[0])
        if error<.0005:break
    return weights,dict(iterations=iteration+1,max_category_fraction_difference=float(error),
        balanced=bool(error<.002),max_weight=float(weights.max()),p99_weight=float(np.percentile(weights,99)),dimensions=len(dims))


def uncertainty(ref, target, weights, seed):
    # Conditional on balancing weights; ordinary multinomial run bootstrap.
    # Same run gets the same multiplier when cohorts overlap in 2011.
    runs,inverse=np.unique(np.r_[P['run'][ref],P['run'][target]],return_inverse=True)
    nr=len(ref);rng=np.random.RandomState(seed);gain=[];deltas=[]
    pairs=[]
    for ix in (ref,target):
        pairs.append([P[key][ix]-P['rich_beta'][ix] for key in ('raw_beta_rich','new_beta_rich')])
    order=[[np.argsort(x,kind='stable') for x in pair] for pair in pairs]
    def width_sorted(x,v,o):
        x=x[o];v=v[o];ok=v>0;x=x[ok];v=v[ok]
        q=np.interp([.16,.84],(np.cumsum(v)-.5*v)/v.sum(),x)
        return (q[1]-q[0])/2
    for _ in range(REPEATS):
        draws=rng.multinomial(len(runs),np.full(len(runs),1/len(runs)))
        rw=draws[inverse[:nr]].astype(float);tw=weights*draws[inverse[nr:]]
        a=[width_sorted(x,rw,o) for x,o in zip(pairs[0],order[0])]
        b=[width_sorted(x,tw,o) for x,o in zip(pairs[1],order[1])]
        ga=100*(1-a[1]/a[0]);gb=100*(1-b[1]/b[0]);gain.append(gb);deltas.append(gb-ga)
    return dict(repeats=REPEATS,cluster='actual run',balancing_weights_fixed=True,
        gain_95ci=np.percentile(gain,[2.5,97.5]),difference_from_reference_95ci=np.percentile(deltas,[2.5,97.5]))


def analyze_year(year):
    ref0=np.flatnonzero(REF);target0=np.flatnonzero((YEAR==year)&P['reference_mask'])
    records=[]
    for stage,dims in zip(STAGES,DIMS):
        ref,target=supported(ref0,target0,dims)
        weights,diagnostics=balance(ref,target,dims)
        a=measure(ref,np.ones(len(ref)));b=measure(target,weights)
        row=dict(stage=stage,reference=a,target=b,balance=diagnostics,
            reference_coverage=len(ref)/len(ref0),target_coverage=len(target)/len(target0),
            gain_difference_pct=b['gain_pct']-a['gain_pct'])
        if stage in ('kinematics','geometry','energy_inputs','timing_inputs'):
            row['bootstrap']=uncertainty(ref,target,weights,20261010+int(year)*10+STAGES.index(stage))
        records.append(row)
    strict=(C['rich_pb']>=.05)&(C['rich_error']<.003)&(C['rich_pmt']>=4)&(C['rich_used']>=6)
    ref,target=supported(ref0[strict[ref0]],target0[strict[target0]],DIMS[1])
    weights,diagnostics=balance(ref,target,DIMS[1])
    strict_result=dict(reference=measure(ref,np.ones(len(ref))),target=measure(target,weights),balance=diagnostics)
    print('Year',year,'raw',round(records[0]['target']['gain_pct'],2),
          'kinematic',round(records[1]['target']['gain_pct'],2),'geometry',round(records[3]['target']['gain_pct'],2),
          'energy',round(records[4]['target']['gain_pct'],2),'timing',round(records[5]['target']['gain_pct'],2),flush=True)
    return dict(year=int(year),n=len(target0),stages=records,strict_reference_sensitivity=strict_result,
        covariate_medians={k:float(np.median(C[k][target0])) for k in ['rich_error','rich_npe_expected','relR','chi_x','chi_y']},
        model_correction_u=dict(mean=float(np.mean(P['raw_u'][target0]-P['corrected_u'][target0])),
        sd=float(np.std(P['raw_u'][target0]-P['corrected_u'][target0]))))


def main():
    global P,C,YEAR,REF,DIMS,REPEATS
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4);parser.add_argument('--repeats',type=int,default=200);args=parser.parse_args();REPEATS=args.repeats
    w.check_freeze();meta=json.loads((OUT/'diagnostic_covariates.json').read_text())
    assert meta['predictions_sha256']==w.ev.digest(OUT/'predictions.npz')
    assert meta['covariates_sha256']==w.ev.digest(OUT/'diagnostic_covariates.npz')
    with np.load(str(OUT/'predictions.npz')) as f:P={k:f[k] for k in f.files}
    with np.load(str(OUT/'diagnostic_covariates.npz')) as f:C={k:f[k] for k in f.files}
    for key in ('run','event','source'):assert np.array_equal(P[key],C[key])
    with np.load(str(w.DEST/'predictions.npz')) as f:oldruns=f['run'][f['reference_mask']]
    REF=P['reference_mask']&(P['run']>=oldruns.min())&(P['run']<=oldruns.max())
    years=np.arange(2011,2027);limits=[datetime(int(y),1,1,tzinfo=timezone.utc).timestamp() for y in years]
    YEAR=years[np.searchsorted(limits,P['run'],side='right')-1]
    DIMS=build_dimensions()
    with mp.Pool(args.workers) as pool:rows=list(pool.imap(analyze_year,np.unique(YEAR)))
    w.ev.dump(OUT/'time_attribution.json',dict(complete=True,beta_model_refitted=False,n=len(P['run']),
        reference_n=int(REF.sum()),reference_run_range=[int(oldruns.min()),int(oldruns.max())],rows=rows,
        controls='Measured-only categorical iterative proportional fitting; no residual/gain used in weights',
        uncertainty='Multinomial actual-run bootstrap, balancing weights fixed',
        limitations='Observational calendar association; energy/timing controls can be mediators. Balancing is marginal plus joint R/beta/radiator strata, not all high-dimensional interactions.',
        covariates_sha256=meta['covariates_sha256'],predictions_sha256=meta['predictions_sha256'],script_sha256=w.ev.digest(Path(__file__))))
    fig,ax=plt.subplots(figsize=(10,6))
    for stage,label,color in [('unadjusted','Observed','#888888'),('kinematics','R, beta, radiator balanced','#5790fc'),
        ('geometry','+ reference quality, geometry','#f89c20'),('energy_inputs','+ energy-deposit inputs','#4caf50'),
        ('timing_inputs','+ timing inputs (mediation sensitivity)','#aa4499')]:
        points=[next(x for x in row['stages'] if x['stage']==stage) for row in rows]
        ax.plot([r['year'] for r in rows],[x['target']['gain_pct'] for x in points],'o-',color=color,label=label)
        if stage=='geometry':
            for row,point in zip(rows,points):ax.plot([row['year']]*2,point['bootstrap']['gain_95ci'],color=color,alpha=.6)
    ax.axhline(0,color='black',linestyle='--',linewidth=1);ax.set_xlabel('Calendar year (UTC)');ax.set_ylabel('Beta residual W68 reduction (%)')
    ax.set_xticks(np.arange(2011,2026,2));ax.legend(frameon=False,fontsize=9);fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(str(OUT/'figures'/('time_attribution.'+ext)),dpi=170)
    plt.close(fig)


if __name__=='__main__':main()
