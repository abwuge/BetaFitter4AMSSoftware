#!/usr/bin/env python3
"""Evaluate the fixed beta model with reference-quality v2, without refitting."""
import argparse
from pathlib import Path
import numpy as np
import rich_transport as rt
import evaluate_rich_transport as ev


def main():
    ap=argparse.ArgumentParser();ap.add_argument('kind',choices=['mc','data']);args=ap.parse_args();kind=args.kind
    dest=rt.ROOT/'results/rich_transport'
    with np.load(str(dest/(kind+'_predictions.npz'))) as f:p={k:f[k] for k in f.files}
    if kind=='mc':
        with np.load(str(dest/'mc_rich_charge.npz')) as f:
            if not np.array_equal(f['run'],p['run']) or not np.array_equal(f['event'],p['event']):raise RuntimeError('Charge identity mismatch')
            q=f['rich_q'][:,1]
    else:
        with np.load('/tmp/betaEstimator-rich-transport/data_evaluation_input.npz') as f:
            if not np.array_equal(f['run'],p['run']) or not np.array_equal(f['event'],p['event']):raise RuntimeError('Flight identity mismatch')
            q=f['rich_q'][:,1]
    loose=p['rich_comparison_sample']
    charge=np.isfinite(q)&(q>4)&(q<8)
    rich=loose&charge
    b0=p['raw_beta_planes'][:,-1];b1=p['new_beta_planes'][:,-1];ref=p['rich_beta']
    label=p['label'] if kind=='mc' else None
    def scored(mask):
        if mask.sum()<2:return {'n':int(mask.sum())}
        result=ev.compare(b0[mask],b1[mask],ref[mask])
        result['bootstrap']=ev.bootstrap(b0[mask],b1[mask],ref[mask],p['source'][mask],label[mask] if label is not None else None)
        return result
    result={'kind':kind,'loose_n':int(loose.sum()),'charge_consistent_n':int(rich.sum()),
            'removed_n':int((loose&~charge).sum()),'rich_all':scored(rich)}
    result['groups']={}
    masks={'NaF':p['NaF'],'aerogel':~p['NaF'],'strict_ring_quality':p['rich_strict'],
           'charge_5_to_7':np.isfinite(q)&(q>5)&(q<7),
           'beta_0.78_0.85':(ref>=.78)&(ref<.85),'beta_0.85_0.90':(ref>=.85)&(ref<.90),
           'beta_0.90_0.95':(ref>=.90)&(ref<.95),'beta_0.95_1.02':(ref>=.95)&(ref<=1.02)}
    if kind=='mc':masks.update({'C12':label==0,'C13':label==1})
    for name,mask in masks.items():result['groups'][name]=scored(rich&mask)
    result['expected_total_loss_GeV']=ev.stats(p['new_loss_GeV'][rich].sum(axis=1))
    result['beta_drop']=ev.stats(1/p['u_new'][rich]-b1[rich])
    result['unpropagated_center_comparison']=ev.compare(1/p['u_raw'][rich],1/p['u_new'][rich],ref[rich])
    saved={'run':p['run'],'event':p['event'],'rich_charge':q,'charge_consistent_reference':rich}
    if kind=='mc':
        tv=p['truth_valid']&np.isfinite(p['u_true'])&(p['u_true']>0)
        true=1/p['u_true'];matched=rich&tv
        result['matched_RICH']=scored(matched)
        result['matched_center_truth']=ev.compare(1/p['u_raw'][matched],1/p['u_new'][matched],true[matched])
        result['matched_center_truth']['bootstrap']=ev.bootstrap(1/p['u_raw'][matched],1/p['u_new'][matched],true[matched],p['source'][matched],label[matched])
        e0=1/p['u_raw'][matched]-true[matched];e1=1/p['u_new'][matched]-true[matched]
        r0=b0[matched]-ref[matched];r1=b1[matched]-ref[matched]
        result['eventwise_agreement']={'raw_error_correlation':float(np.corrcoef(e0,r0)[0,1]),
                                      'new_error_correlation':float(np.corrcoef(e1,r1)[0,1]),
                                      'fraction_improved_by_truth':float(np.mean(np.abs(e1)<np.abs(e0))),
                                      'fraction_improved_by_RICH':float(np.mean(np.abs(r1)<np.abs(r0))),
                                      'improvement_decision_agreement':float(np.mean((np.abs(e1)<np.abs(e0))==(np.abs(r1)<np.abs(r0))))}
        with np.load('/tmp/betaEstimator-rich-transport/mc_evaluation_input.npz') as f:
            mom=f['eval_mevmom1'].astype(float);xyz=f['eval_mevcoo1']
            trmom=f['eval_mtrmom'];pri=f['eval_mtrpri'];pos=f['tk_pos'];direction=f['tk_dir']
        mass=np.where(label==0,11.17793,12.11221);energy17=np.sqrt(mom[:,17]**2+mass*mass);bp=mom[:,17]/energy17
        ok=(mom[:,17]>0)&np.isfinite(mom[:,17])&(np.abs(xyz[:,17,2]+69.975)<.06)&(pri[:,4]==1)&(trmom[:,4]>0)
        result['upstream_truth_checks']={}
        for name,mask in [('loose',loose),('charge_consistent',rich),('charge_consistent_primary_L9',rich&(pri[:,8]==1)&(trmom[:,8]>0)&(trmom[:,8]<=mom[:,17]+.01))]:
            sel=mask&tv&ok;delta=ref[sel]-bp[sel]
            result['upstream_truth_checks'][name]={'n':int(sel.sum()),'above_upstream_truth_gt_0.005':int((delta>.005).sum()),
                                                   'above_upstream_truth_gt_0.02':int((delta>.02).sum()),
                                                   'below_upstream_truth_lt_minus0.02':int((delta<-.02).sum())}
        sel=matched&ok;ix=np.flatnonzero(sel)
        material=rt.MaterialMap();geo=material.build_columns(pos[ix],direction[ix]);del material
        exact=rt.propagate(p['u_true'][ix],mom[ix,10]/6,geo)
        valid=exact['valid'];ix=ix[valid]
        exact_rich=exact['beta_planes'][valid,-1];exact_post=exact['beta_planes'][valid,6]
        result['exact_center_transport_diagnostic']={
            'n':len(ix),'RICH_prediction_minus_measured_beta':ev.stats(exact_rich-ref[ix]),
            'postS4_prediction_minus_truth_beta':ev.stats(exact_post-bp[ix]),
            'fraction_RICH_abs_error_gt_0.005':float(np.mean(np.abs(exact_rich-ref[ix])>.005)),
            'fraction_RICH_abs_error_gt_0.02':float(np.mean(np.abs(exact_rich-ref[ix])>.02))}
        saved['exact_center_rich_beta']=np.full(len(q),np.nan);saved['exact_center_rich_beta'][ix]=exact_rich
        saved['true_beta_postS4']=np.where(sel,bp,np.nan)
        result['transport_sensitivity_matched_clean']={}
        baseidx=np.flatnonzero(sel)
        for name,kw in [('material_0.8',{'material_scale':.8}),('material_1.2',{'material_scale':1.2}),
                        ('mass_C12',{'mass_override':11.17793}),('mass_C13',{'mass_override':12.11221})]:
            v0=rt.propagate(p['u_raw'][baseidx],p['R'][baseidx],geo,**kw)
            v1=rt.propagate(p['u_new'][baseidx],p['R'][baseidx],geo,**kw)
            m=v0['valid']&v1['valid']
            result['transport_sensitivity_matched_clean'][name]=ev.compare(v0['beta_planes'][m,-1],v1['beta_planes'][m,-1],ref[baseidx][m])
    np.savez_compressed(str(dest/(kind+'_reference_v2_masks.npz')),**saved)
    ev.dump(dest/(kind+'_reference_v2_results.json'),result)
    print(kind,'RICH v2',int(rich.sum()),'/',int(loose.sum()),flush=True)
    keys=['rich_all','matched_RICH','matched_center_truth','eventwise_agreement','upstream_truth_checks','exact_center_transport_diagnostic']
    print(__import__('json').dumps({k:result[k] for k in keys if k in result},indent=2,default=lambda x:x.tolist()),flush=True)

if __name__=='__main__':main()
