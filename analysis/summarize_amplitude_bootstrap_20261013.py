#!/usr/bin/env python3
"""Paired train-run plus evaluation-run uncertainty for amplitude-noise ablation."""
import argparse,json
from pathlib import Path
import numpy as np
from bootstrap_amplitude_noise_20261013 import BOOT,OUT
from research_conditional_mc import OUT as MCOUT,stamp
from research_flight_conditional import load,write,digest
from evaluate_amplitude_noise_20261013 import join,identity
from evaluate_conditional_mc import metrics,pair_delta

SCORE=OUT/'scoring'
LABELS=[family+'_'+variant for family in ['native_gaussianR','native_studentR5'] for variant in ['geometry','readout']]


def prepare():
    path=SCORE/'joint_bootstrap_protocol.json'
    if path.exists():raise FileExistsError('Bootstrap evaluation protocol exists')
    write(path,dict(created_utc=stamp(),training_replicates=128,evaluation_draws_per_training_replicate=16,seed=195202610,
        training='Same paired whole-training-run draw for all four models; per-draw feature normalization and two complete optimization starts',
        evaluation='Sixteen independent whole-evaluation-run draws per frozen resampled calibration; same event draws across candidate comparisons',
        interpretation='128 distinct recalibrations and 2048 nested train/evaluation combinations, NOT 2048 independent training replicates',
        spaces=['u','beta'],quantiles=[2.5,50,97.5],comparisons='Four models versus native, readout versus same-family same-resample geometry',
        bootstrap_protocol_sha256=digest(BOOT/'protocol.json'),reference_sha256=digest(MCOUT/'scoring/reference.npz'),script_sha256=digest(Path(__file__)),
        failures='Include all completed repetitions; convergence failures prevent acceptance and must be reported, never removed by beta error',RICH_read=False))
    print('AMPLITUDE_JOINT_BOOT_PROTOCOL',flush=True)


def score():
    protocol=json.loads((SCORE/'joint_bootstrap_protocol.json').read_text());boot=json.loads((BOOT/'protocol.json').read_text())
    if digest(Path(__file__))!=protocol['script_sha256'] or digest(BOOT/'protocol.json')!=protocol['bootstrap_protocol_sha256']:raise ValueError('Bootstrap evaluation inputs changed')
    if (SCORE/'joint_bootstrap.json').exists():raise FileExistsError('Joint bootstrap already summarized')
    notes=[];failures=[];n=protocol['training_replicates']
    for rep in range(n):
        path=BOOT/('%04d.json'%rep);note=json.loads(path.read_text())
        if not note['complete'] or note['replicate']!=rep or note['protocol_sha256']!=protocol['bootstrap_protocol_sha256'] or digest(path.with_suffix('.npz'))!=note['output_sha256']:raise ValueError('Invalid completed bootstrap '+str(rep))
        if sorted(r['label'] for r in note['models'])!=sorted(LABELS):raise ValueError('Missing bootstrap model')
        for row in note['models']:
            if not row['model']['success'] or row['model']['projected_gradient_max']>boot['nll_convergence_tolerance']:
                failures.append(dict(replicate=rep,label=row['label'],success=row['model']['success'],pg=row['model']['projected_gradient_max']))
        notes.append(note)
    if digest(MCOUT/'scoring/reference.npz')!=protocol['reference_sha256']:raise ValueError('Bootstrap reference changed')
    first=load(BOOT/'0000.npz');ref=load(MCOUT/'scoring/reference.npz');a=load(MCOUT/'measured.npz');ir=join(first,ref);im=join(first,a);valid=ref['valid'][ir]
    truth=ref['u'][ir][valid];native=1/a['tof_betah'][im][valid];runs=first['run'][valid];groups=[np.flatnonzero(runs==r) for r in np.unique(runs)]
    pairs=[(label,'native') for label in LABELS]+[(family+'_readout',family+'_geometry') for family in ['native_gaussianR','native_studentR5']]
    store={};raw={}
    for space in protocol['spaces']:
        for label,old in pairs:
            key=space+'__'+label+'__vs__'+old;store[key]=dict(train_only=[],joint=[])
    original={label:load(OUT/'mc'/label/'predictions_before_truth.npz')['u'][valid] for label in LABELS};original['native']=native
    for rep,note in enumerate(notes):
        data=load(BOOT/('%04d.npz'%rep))
        if not np.array_equal(identity(first),identity(data)):raise ValueError('Bootstrap evaluation identities changed')
        values={label:data[label+'_u'][valid] for label in LABELS};values['native']=native
        rng=np.random.RandomState(protocol['seed']+rep);draws=[np.concatenate([groups[j] for j in rng.randint(0,len(groups),len(groups))]) for _ in range(protocol['evaluation_draws_per_training_replicate'])]
        for space in protocol['spaces']:
            err={k:v-truth if space=='u' else 1/v-1/truth for k,v in values.items()}
            for label,old in pairs:
                key=space+'__'+label+'__vs__'+old
                store[key]['train_only'].append(pair_delta(err[old],err[label]))
                store[key]['joint'].extend(pair_delta(err[old][ix],err[label][ix]) for ix in draws)
        if (rep+1)%32==0:print('AMPLITUDE_BOOT_SCORING',rep+1,'/',n,flush=True)
    summary={}
    for space in protocol['spaces']:
        err={k:v-truth if space=='u' else 1/v-1/truth for k,v in original.items()};summary[space]={}
        for label,old in pairs:
            key=space+'__'+label+'__vs__'+old;rows=store[key]
            out=dict(point=pair_delta(err[old],err[label]),train_only={},joint={})
            for kind,rep in rows.items():
                for metric in rep[0]:
                    values=np.array([r[metric] for r in rep]);raw[key+'__'+kind+'__'+metric]=values;out[kind][metric]=np.percentile(values,protocol['quantiles'])
            summary[space][label+'_versus_'+old]=out
    np.savez_compressed(SCORE/'joint_bootstrap_metrics.npz',**raw)
    result=dict(complete=True,created_utc=stamp(),events=len(truth),runs=len(groups),replicates=n,nested_combinations=n*protocol['evaluation_draws_per_training_replicate'],
        convergence_failures=failures,accepted_numerical_runs=not failures,intervals=summary,protocol_sha256=digest(SCORE/'joint_bootstrap_protocol.json'),
        metrics_sha256=digest(SCORE/'joint_bootstrap_metrics.npz'),reference_sha256=protocol['reference_sha256'],
        evidence_sha256={str((BOOT/('%04d.json'%i)).relative_to(OUT)):digest(BOOT/('%04d.json'%i)) for i in range(n)},RICH_read=False)
    write(SCORE/'joint_bootstrap.json',result)
    for space in protocol['spaces']:
        for label,old in pairs:
            row=summary[space][label+'_versus_'+old]
            print('AMPLITUDE_BOOT_INTERVAL',space,label,old,'width',row['joint']['w68_improvement'],'rmse',row['joint']['rmse_improvement'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','score']);args=parser.parse_args()
    if args.stage=='prepare':prepare()
    else:score()
