#!/usr/bin/env python3
"""Additional source/run holdout, measurement preparation separate from truth."""
import argparse,hashlib,json,subprocess,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import uproot
import research_conditional_mc as previous
import research_flight_material as material
import feature_conditional as feature
import conditional_readout_population as population
import joint_adc_likelihood as joint
import free_timing_transport as transport
from research_amplitude_noise_20261013 import OUT as AMPOUT,feature_arrays
from research_readout_population import OUT as POPOUT
from research_joint_adc import OUT as JOINTOUT
from research_flight_conditional import load,write,digest
from research_conditional_mc_reference import RAW_C12,ENV
from research_conditional_mc import ROOT,OUT as MCOUT,stamp

OUT=ROOT/'results/flight_calibration_research/continuation_20261013/additional_mc'
SOURCE=Path('/lustre02/user/yachen/public/hxwu/amsd69n_mc_z1_z2_z6_z8_20260903/Z6_C12')


def plan():
    OUT.mkdir(exist_ok=True)
    if (OUT/'protocol.json').exists():raise FileExistsError('Additional-run plan already frozen')
    prior=json.loads((MCOUT/'protocol.json').read_text());used={r['path'] for r in prior['records']}
    available=[p for p in SOURCE.glob('*.root') if str(p) not in used]
    available.sort(key=lambda p:hashlib.sha256(('additional_mc_20261013:'+p.name).encode()).digest())
    selected=available[:128]
    records=[]
    for i,path in enumerate(selected):
        raw=RAW_C12/(path.stem.replace('_','.')+'.root')
        if not raw.exists():raise FileNotFoundError('Prespecified source lacks native reference: '+str(raw))
        records.append(dict(source=i,path=str(path),raw=str(raw),entry_stop=8000))
    write(OUT/'protocol.json',dict(created_utc=stamp(),records=records,available_additional_sources=len(available),maximum_events=20000,
        source_selection='First 128 files by SHA256(additional_mc_20261013:filename), outside the complete prior conditional-MC input manifest; no reference or reconstructed-error choice',
        event_selection='Exact existing measured carbon/native-four-used/Kalman-quality domain plus positive end ADC; remove every run in the prior measured MC cache before inference',
        sample_downselection='Deterministic SHA256(additional_event:run:event) to at most 20000 measured-selected events',
        model_domain='All previous amplitude covariance, conditional-population and new joint-ADC models, frozen unchanged; no normalization or calibration refit on these additional events',
        reference='Same previously audited native parent decoder and primary plane10 L5 conditions, accessed only after all output arrays are frozen',
        reference_conditions=dict(plane=10,z=1.698,z_tolerance=.05,track_distance=1.,charge=6,downward=True),
        evaluation_primary='native_studentR5_joint_transport_u versus stored native beta; other methods are prespecified controls',
        evaluation_spaces=['u','beta'],evaluation_bootstrap=2048,bins=dict(R=[2.5,3,3.5,4.5,6],beta=[0,.8,.85,.9,.95,1.01]),
        development_caveat='Not used in this conditional-model fitting/design data manifest; entire historical project exposure is not certified, so not called universally unseen or a new flight blind test',
        previous_protocol_sha256=digest(MCOUT/'protocol.json'),previous_measurement_sha256=digest(MCOUT/'measured.npz'),
        measurement_selector_sha256=digest(ROOT/'analysis/research_conditional_mc.py'),script_sha256=digest(Path(__file__)),
        MC_truth_read=False,RICH_read=False))
    print('ADDITIONAL_MC_PLAN',len(records),'of',len(available),flush=True)


def prepare():
    p=json.loads((OUT/'protocol.json').read_text())
    if (OUT/'measured.npz').exists():raise FileExistsError('Additional measurements already saved')
    if digest(Path(__file__))!=p['script_sha256'] or digest(ROOT/'analysis/research_conditional_mc.py')!=p['measurement_selector_sha256']:raise ValueError('Additional measurement preparation changed')
    old=load(MCOUT/'measured.npz');oldruns=np.unique(old['run']);pieces=[];audits=[]
    def one(record):
        a,note=previous.measured_record(record)
        with uproot.open(record['path']) as stream:
            raw=stream['amstreea'].arrays(['run','event','itrtrack','ibetah','tof_Aadc'],entry_stop=record['entry_stop'],library='np')
        index=a['entry']
        for k in ['run','event','itrtrack','ibetah']:
            if not np.array_equal(a[k],raw[k][index]):raise ValueError('Additional ADC identity mismatch')
        adc=raw['tof_Aadc'][index].astype(float);valid=np.isfinite(adc).all((1,2))&(adc>0).all((1,2))
        overlap=np.isin(a['run'],oldruns);keep=valid&~overlap
        a['adc']=adc;a={k:v[keep] for k,v in a.items()}
        note.update(ADC_branch='tof_Aadc',removed_ADC_invalid=int((~valid).sum()),removed_previous_run_overlap=int(overlap.sum()),common_domain_events=int(keep.sum()),truth_read=False)
        return a,note
    with ThreadPoolExecutor(max_workers=4) as pool:
        for i,(a,note) in enumerate(pool.map(one,p['records'])):
            pieces.append(a);audits.append(note)
            if (i+1)%16==0:print('ADDITIONAL_MC_READ',i+1,'/',len(p['records']),flush=True)
    a={k:np.concatenate([v[k] for v in pieces]) for k in pieces[0]}
    ident=material.keys(a)
    if len(np.unique(ident))!=len(ident):raise ValueError('Additional duplicate run/event')
    count=len(ident);rank=np.array([int.from_bytes(hashlib.sha256(('additional_event:%d:%d'%(r,e)).encode()).digest()[:8],'big') for r,e in zip(a['run'],a['event'])],np.uint64)
    selected=np.sort(np.argsort(rank)[:min(p['maximum_events'],count)]);a={k:v[selected] for k,v in a.items()};a['fold']=np.zeros(len(selected),np.int8)
    if np.intersect1d(a['run'],oldruns).size:raise ValueError('Additional-run separation failed')
    np.savez_compressed(OUT/'measured.npz',**a)
    values=np.column_stack([a['run'],a['event'],a['R'],a['bar'],a['x']*previous.C,a['tk_pos'].reshape(len(selected),-1),a['tk_dir'].reshape(len(selected),-1)])
    np.savetxt(OUT/'geometry_input.txt',values,fmt=['%d','%d','%.12g']+['%d']*4+['%.12g']*58)
    write(OUT/'preparation.json',dict(complete=True,created_utc=stamp(),events=len(selected),parent_selected=count,runs=len(np.unique(a['run'])),
        prior_run_overlap=0,records=audits,measured_sha256=digest(OUT/'measured.npz'),geometry_input_sha256=digest(OUT/'geometry_input.txt'),RICH_read=False,MC_truth_read=False))
    print('ADDITIONAL_MC_READY',len(selected),len(np.unique(a['run'])),flush=True)


def geometry():
    if (OUT/'geometry.root').exists() or (OUT/'moments.npz').exists():raise FileExistsError('Additional native geometry already exists')
    executable=ROOT/'results/flight_calibration_research/mean_loss_likelihood/supplement'
    command=['bash','-c','source "$1" >/dev/null; exec "$2" "$3" "$4"','geometry',str(ENV),str(executable),str(OUT/'geometry_input.txt'),str(OUT/'geometry.root')]
    start=time.monotonic()
    with (OUT/'geometry.log').open('w') as stream:r=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
    if r.returncode:raise RuntimeError('Additional native geometry failed')
    original=material.OUT;material.OUT=OUT
    try:material.moments()
    finally:material.OUT=original
    g=load(OUT/'moments.npz')
    if not g['valid'].all():raise ValueError('Additional geometry invalid; do not outcome-filter')
    write(OUT/'geometry_execution.json',dict(complete=True,created_utc=stamp(),elapsed_seconds=time.monotonic()-start,executable_sha256=digest(executable),
        source_sha256=digest(ROOT/'analysis/flight_material_geometry.C'),geometry_sha256=digest(OUT/'geometry.root'),moments_sha256=digest(OUT/'moments.npz'),RICH_read=False,MC_truth_read=False))


def predict():
    if (OUT/'output_freeze.json').exists():raise FileExistsError('Additional output already frozen')
    a=load(OUT/'measured.npz');g=load(OUT/'moments.npz');choose=np.ones(len(a['run']),bool);adc=a['adc']
    if not np.array_equal(material.keys(a),material.keys(g)) or not g['valid'].all():raise ValueError('Additional prediction geometry mismatch')
    specs=[]
    for family in ['native_gaussianR','native_studentR5']:
        for variant in ['geometry','readout']:specs.append(('amplitude',family+'_'+variant,AMPOUT/'mc'/(family+'_'+variant)))
        for variant in ['population_only','joint']:specs.append(('population',family+'_'+variant,POPOUT/'mc'/(family+'_'+variant)))
        specs.append(('joint_adc',family,JOINTOUT/'mc'/family))
    frozen={}
    for kind,label,directory in specs:
        meta=json.loads((directory/'measurement_evaluation.json').read_text())
        if not meta['complete'] or digest(directory/'freeze.json')!=meta['freeze_sha256'] or digest(directory/'predictions_before_truth.npz')!=meta['output_sha256']:raise ValueError('Additional inference candidate not frozen')
        frozen[kind+'/'+label]=dict(model=meta['freeze_sha256'],original_output=meta['output_sha256'],model_path=str(directory/'freeze.json'))
    output={k:a[k] for k in ['run','event','source','entry','R']};xc=g['x']-g['x'].mean(1)[:,None];w0=xc/np.sum(xc*xc,axis=1)[:,None]
    output['ols_u']=np.sum(w0*a['t_native'],axis=1);output['native_u']=1/a['tof_betah'];output['constant_u']=np.ones(len(choose))
    ap=json.loads((AMPOUT/'mc/protocol.json').read_text());pp=json.loads((POPOUT/'mc/protocol.json').read_text());jp=json.loads((JOINTOUT/'mc/protocol.json').read_text())
    for kind,label,directory in specs:
        theta=np.asarray(json.loads((directory/'freeze.json').read_text())['model']['theta'])
        if kind=='amplitude':
            variant='geometry' if label.endswith('_geometry') else 'readout';ff,_=feature_arrays(a,adc,choose,choose,variant,ap['variants'][variant]['normalization'])
            w,Q=feature.precision_weights(theta,g['x'],ff);output['amplitude_'+label+'_u']=np.sum(w*a['t_native'],axis=1)
        elif kind=='population':
            ff,_=feature_arrays(a,adc,choose,choose,'readout',pp['feature_normalization']);w,Q=population.precision_weights(theta,g['x'],ff)
            output['population_'+label+'_u']=np.sum(w*a['t_native'],axis=1)
        else:
            ff,_=joint.response_features(a,choose,jp['feature_normalization']);yy=joint.observation(a,adc);result=joint.reconstruct(theta,g['x'],yy,ff)
            for branch in ['time_u','adc_u','joint_u']:output[label+'_'+branch]=result[branch]
            for branch,weight in [('time',result['time_weights']),('joint',result['weights'][:,:3].dot(previous.f.D))]:
                solved=transport.solve_center(result[branch+'_u'],a['R'],weight,g['moment'],theta[:2])
                if not solved['converged'].all():raise ValueError('Additional mean-loss inverse failed')
                output[label+'_'+branch+'_transport_u']=solved['u']
            output[label+'_time_injection_gain']=result['time_injection_gain']
    if not all(np.isfinite(v).all() for k,v in output.items() if k.endswith('_u')):raise ValueError('Nonfinite additional output; no filtering allowed')
    np.savez_compressed(OUT/'predictions_before_truth.npz',**output)
    write(OUT/'output_freeze.json',dict(complete=True,created_utc=stamp(),events=len(choose),models=frozen,branches=[k for k in output if k.endswith('_u')],
        output_sha256=digest(OUT/'predictions_before_truth.npz'),measured_sha256=digest(OUT/'measured.npz'),moments_sha256=digest(OUT/'moments.npz'),
        protocol_sha256=digest(OUT/'protocol.json'),normalization_refit=False,calibration_refit=False,RICH_read=False,MC_truth_read=False))
    print('ADDITIONAL_OUTPUT_FROZEN',len(choose),len(frozen),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['plan','prepare','geometry','predict']);args=p.parse_args()
    globals()[args.stage]()
