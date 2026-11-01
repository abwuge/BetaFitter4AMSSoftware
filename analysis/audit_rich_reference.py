#!/usr/bin/env python3
"""Post-freeze reference audit and exactly matched MC evaluation populations."""
import json
from pathlib import Path
import numpy as np
import rich_transport as rt
import evaluate_rich_transport as ev


def main():
    dest=rt.ROOT/'results/rich_transport'
    with np.load(str(dest/'mc_predictions.npz')) as f:p={k:f[k] for k in f.files}
    with np.load('/tmp/betaEstimator-rich-transport/mc_evaluation_input.npz') as f:a={k:f[k] for k in f.files}
    r=p['rich_comparison_sample'];tv=p['truth_valid']&np.isfinite(p['u_true'])&(p['u_true']>0)
    matched=r&tv
    ref=p['rich_beta'];true=1/p['u_true']
    b0=p['raw_beta_planes'][:,-1];b1=p['new_beta_planes'][:,-1]
    output={'matched_center_truth':ev.compare(1/p['u_raw'][matched],1/p['u_new'][matched],true[matched]),
            'matched_RICH':ev.compare(b0[matched],b1[matched],ref[matched]),
            'matched_RICH_bootstrap':ev.bootstrap(b0[matched],b1[matched],ref[matched],p['source'][matched],p['label'][matched])}
    mass=np.where(p['label']==0,11.17793,12.11221)
    mom=a['eval_mevmom1'].astype(float);xyz=a['eval_mevcoo1'].astype(float)
    ep=np.sqrt(mom[:,17]**2+mass*mass)
    bp=mom[:,17]/ep
    post=(mom[:,17]>0)&np.isfinite(mom[:,17])&(np.abs(xyz[:,17,2]+69.975)<.06)
    primary5=(a['eval_mtrpri'][:,4]==1)&(a['eval_mtrmom'][:,4]>0)
    primary9=(a['eval_mtrpri'][:,8]==1)&(a['eval_mtrmom'][:,8]>0)
    seq=(a['eval_mtrmom'][:,8] <= mom[:,17]+.01)
    masks={'all':np.ones(len(r),bool),'NaF':p['NaF'],'aerogel':~p['NaF'],
           'C12':p['label']==0,'C13':p['label']==1,'strict':p['rich_strict'],
           'primary_through_L9':primary9&seq,'NaF_primary_through_L9':p['NaF']&primary9&seq}
    output['upstream_truth_reference_checks']={}
    for name,mask in masks.items():
        s=r&tv&post&primary5&mask
        if s.sum()<2:
            output['upstream_truth_reference_checks'][name]={'n':int(s.sum())};continue
        diff=ref[s]-bp[s]
        rec={'n':int(s.sum()),'RICH_minus_upstream_true_beta':ev.stats(diff),
             'RICH_above_upstream_truth_gt_0.005':int(np.sum(diff>.005)),
             'RICH_above_upstream_truth_gt_0.02':int(np.sum(diff>.02)),
             'RICH_below_upstream_truth_lt_minus0.02':int(np.sum(diff<-.02)),
             'RICH_resolution_comparison':ev.compare(b0[s],b1[s],ref[s]),
             'center_truth_comparison':ev.compare(1/p['u_raw'][s],1/p['u_new'][s],true[s])}
        for k in range(3):
            rec['stored_RICH_beta%d_minus_postS4_truth'%k]=ev.stats(a['rich_beta'][s,k]-bp[s])
        output['upstream_truth_reference_checks'][name]=rec
    impossible=r&tv&post&primary5&primary9&seq&(ref>bp+.02)
    indices=np.flatnonzero(impossible)
    # Store event identities and all RICH beta variants for the strongest
    # counterexamples, without changing the main validation selection.
    order=indices[np.argsort((ref-bp)[indices])[::-1]][:12]
    output['primary_survives_examples']=[{
        'run':int(p['run'][i]),'event':int(p['event'][i]),'label':int(p['label'][i]),'source':int(p['source'][i]),
        'NaF':bool(p['NaF'][i]),'R_GV':float(p['R'][i]),'beta_true_center':float(true[i]),
        'beta_true_postS4':float(bp[i]),'beta_rich_variants':a['rich_beta'][i].tolist(),
        'p_true_postS4_GeV':float(mom[i,17]),'p_primary_L9_GeV':float(a['eval_mtrmom'][i,8]),
        'rich_prob':float(a['rich_pb'][i]),'rich_used':int(a['rich_used'][i]),
        'rich_pmt':int(a['rich_pmt'][i]),'rich_beta_error':float(a['rich_ebeta1'][i]),
        'rich_beta_consistency':float(a['rich_BetaConsistency'][i]),
        'rich_pmt_charge_consistency':float(a['rich_PMTChargeConsistency'][i])} for i in order]
    ev.dump(dest/'mc_reference_audit.json',output)
    summary={'matched_center_truth':output['matched_center_truth'],'matched_RICH':output['matched_RICH'],
             'matched_RICH_bootstrap':output['matched_RICH_bootstrap'],
             'reference_flags':{k:{q:v[q] for q in ['n','RICH_above_upstream_truth_gt_0.005','RICH_above_upstream_truth_gt_0.02','RICH_below_upstream_truth_lt_minus0.02'] if q in v}
                                for k,v in output['upstream_truth_reference_checks'].items()},
             'examples':output['primary_survives_examples'][:2]}
    print(json.dumps(summary,indent=2,default=lambda x:x.tolist()),flush=True)

if __name__=='__main__':main()
