#!/usr/bin/env python3
"""Score the already frozen curve sensitivity on the old MC development pool."""
import json
from pathlib import Path
import numpy as np
from audit_free_curve_timing import DEST
from research_conditional_mc import OUT as MC,stamp
from research_flight_conditional import load,write,digest
from evaluate_additional_mc_20261013 import matrix_metrics,matrix_difference,METRICS,DELTA


def main():
    target=DEST/'mc_scoring_protocol.json'
    if target.exists():raise FileExistsError('Curve MC scoring was already started')
    audit=json.loads((DEST/'evaluation.json').read_text())
    if not audit['complete'] or not audit['no_beta_reference_read']:raise ValueError('Complete the measured-only curve audit first')
    families=['native_gaussianR','native_studentR5'];paths=[DEST/('mc_'+name+'.npz') for name in families]
    reference_path=MC/'scoring/reference.npz'
    write(target,dict(created_utc=stamp(),audit_sha256=digest(DEST/'evaluation.json'),source_sha256=digest(Path(__file__)),
        predictions={str(path):digest(path) for path in paths},reference_sha256=digest(reference_path),
        purpose='Development sensitivity comparison of one fixed numerical refinement, not new training or a blind confirmation.',
        comparisons='Both MC noise-response families, refined output versus their own original continuous material inverse.',
        bootstrap_replicates=4096,seed=19131026,unit='entire paired evaluation run',
        deterministic_bounds='For the same events and any reference, absolute RMSE change <= RMS of output difference, and absolute W68 change <= maximum absolute output difference.',
        no_reference_residual_cut=True,no_candidate_pruning=True))
    reference=load(reference_path);lookup={(int(run),int(event)):i for i,(run,event) in enumerate(zip(reference['run'],reference['event']))}
    if len(lookup)!=len(reference['run']):raise ValueError('Duplicate reference identity')
    records=[]
    for family,path in zip(families,paths):
        pred=load(path);ir=np.asarray([lookup[(int(run),int(event))] for run,event in zip(pred['run'],pred['event'])])
        valid=reference['valid'][ir];truth=reference['u'][ir][valid];runs=pred['run'][valid]
        if not pred['converged'][valid].all():raise ValueError('Report nonconverged curve predictions rather than removing them')
        values=np.column_stack([pred['original_u'][valid],pred['u'][valid]])
        if not np.isfinite(values).all() or np.any(values==0):raise ValueError('Invalid output on fixed cohort')
        errors={'u':values-truth[:,None],'beta':1/values-1/truth[:,None]}
        groups=[np.flatnonzero(runs==run) for run in np.unique(runs)];rng=np.random.RandomState(19131026)
        bootstrap={space:np.empty((4096,len(DELTA))) for space in errors}
        for rep in range(4096):
            selected=np.concatenate([groups[j] for j in rng.randint(0,len(groups),len(groups))])
            for space,error in errors.items():
                mm=matrix_metrics(error[selected]);bootstrap[space][rep]=matrix_difference(mm[0],mm[1])
        row=dict(family=family,events=len(runs),runs=len(groups),spaces={})
        for space,error in errors.items():
            metrics=matrix_metrics(error);difference=matrix_difference(metrics[0],metrics[1]);change=error[:,1]-error[:,0]
            rms=float(np.sqrt(np.mean(change**2)));maximum=float(np.max(abs(change)))
            if abs(metrics[1,3]-metrics[0,3])>rms+1e-12 or abs(metrics[1,2]-metrics[0,2])>maximum+1e-12:
                raise ValueError('Output perturbation bound failed')
            row['spaces'][space]=dict(original=dict(zip(METRICS,metrics[0])),refined=dict(zip(METRICS,metrics[1])),
                difference=dict(zip(DELTA,difference)),paired_run_interval={key:np.percentile(bootstrap[space][:,k],[2.5,50,97.5]) for k,key in enumerate(DELTA)},
                absolute_RMSE_change_bound=rms,absolute_W68_change_bound=maximum,
                maximum_possible_relative_RMSE_gain=rms/metrics[0,3],maximum_possible_relative_W68_gain=maximum/metrics[0,2])
        records.append(row);print('FREE_CURVE_MC_BETA',family,json.dumps(row['spaces']['beta'],default=lambda a:a.tolist()),flush=True)
    write(DEST/'mc_evaluation.json',dict(complete=True,created_utc=stamp(),rows=records,protocol_sha256=digest(target),
        MC_truth_used_only_for_this_development_evaluation=True,flight_parameters_changed=False,RICH_read=False))


if __name__=='__main__':main()
