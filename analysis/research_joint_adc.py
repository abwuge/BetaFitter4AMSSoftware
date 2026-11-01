#!/usr/bin/env python3
"""Calibrate a joint time / raw end-ADC response independently on MC and flight."""
import argparse,json
from pathlib import Path
import numpy as np
import joint_adc_likelihood as model
import flight_conditional_likelihood as base
import free_timing_transport as transport
from research_amplitude_noise_20261013 import OUT as AMPOUT,inputs,CONFIGS
from research_feature_noise import solve
from research_flight_conditional import load,write,digest
from research_conditional_mc import stamp
ROOT=model.ROOT;OUT=model.OUT


def prepare(dataset):
    directory=OUT/dataset;directory.mkdir(exist_ok=True)
    if (directory/'protocol.json').exists():raise FileExistsError('Joint protocol already frozen')
    folder,a,g,adc,domain,train,test=inputs(dataset);features,norm=model.response_features(a,train)
    seeds={}
    for name,field,nu in CONFIGS[dataset]:
        p=(folder if dataset=='mc' else folder/'conditional')/name/'freeze.json';seeds[name]=dict(path=str(p),sha256=digest(p))
    write(directory/'protocol.json',dict(created_utc=stamp(),dataset=dataset,configs=CONFIGS[dataset],nfeatures=features.shape[2],nparameters=39+4*features.shape[2],
        feature_normalization=norm,features=['intercept']+model.GEOMETRY_NAMES+['bar_%d'%b for b in range(10)],
        geometry_note='Measured Tracker L5 transverse coordinates, not claimed to be a native TOF plane intersection; independent layer coefficients and native associated bar categories',
        observation='Three orthonormal Dt plus four log geometric means of direct associated end ADC; positive-ADC common domain inherited unchanged',
        ADC_mean='layer intercept + measured-geometry/bar response + exp(log_lambda_layer)*(u_L5 - 1.15)',
        ADC_response_scope='Empirical monotone logADC response linear in inverse speed over the selected range; not a calibrated visible-MeV or quenching/Bethe response',
        density='Seven-dimensional Student-5 with full common positive-definite scale, integrated over C12/C13 and true inverse rigidity; normalized over Dt and logADC',
        qrange=[.02,.7],nodes=128,evaluation_nodes=256,bounds=model.bounds(features.shape[2]),penalty=.01,
        penalty_formula='0.5 * 0.01 * sum(non-intercept ADC-response coefficients ** 2); no precision, material, slope or source penalty',
        initializations='Own-dataset timing model with lambda=1 plus independent physical timing initialization with lambda=3; each ADC intercept is own training median logADC',
        optimization='Two starts; up to three consecutive 900-iteration solves per start if needed; select minimum penalized training joint likelihood only',
        event_outputs=['time-only marginal GLS','ADC-only marginal GLS','joint continuous GLS','time-only plus continuous-mass L5 loss inverse','joint plus continuous-mass L5 loss inverse'],
        posterior_boundary='Event outputs do not evaluate an isotope, rigidity or source-spectrum posterior. Mean-loss inverse uses measured rigidity and continuous mass.',
        coherent_response='Joint input shift by b*(Dx,lambda) has exact unit gain; fixed-ADC time-only shift generally has gain different from one and is explicitly reported',
        not_pure_time_estimator=True,selected_response_model_only=True,
        assumptions_not_validated=['selected native-hit/charge truncation','linear logADC-versus-u response','constant seven-dimensional noise scale','full-domain first-order mean loss','period dependence and absolute ADC calibration'],
        evaluation='Same developed MC event pool; no new blind flight RICH scoring; joint 7-D NLL must not be numerically compared with previous 3-D conditional NLL',
        train_events=int(train.sum()),test_events=int(test.sum()),train_runs=np.unique(a['run'][train]),test_runs=np.unique(a['run'][test]),
        MC_truth_read=False,RICH_read=False,MC_parameters_transferred_to_flight=False,seeds=seeds,
        measurement_sha256=digest(folder/'measured.npz'),geometry_sha256=digest(folder/'moments.npz'),amplitude_protocol_sha256=digest(AMPOUT/dataset/'protocol.json'),
        implementation_sha256={n:digest(ROOT/'analysis'/n) for n in ['joint_adc_likelihood.cpp','joint_adc_likelihood.py','research_joint_adc.py','free_timing_transport.py']},
        library_sha256=digest(OUT/'libjointadc.so')))
    print('JOINT_ADC_PROTOCOL',dataset,int(train.sum()),int(test.sum()),flush=True)


def fit(dataset,name):
    folder,a,g,adc,domain,train,test=inputs(dataset);directory=OUT/dataset;p=json.loads((directory/'protocol.json').read_text())
    for f,h in p['implementation_sha256'].items():
        if digest(ROOT/'analysis'/f)!=h:raise ValueError('Joint source changed '+f)
    if digest(OUT/'libjointadc.so')!=p['library_sha256']:raise ValueError('Joint library changed')
    if digest(folder/'measured.npz')!=p['measurement_sha256'] or digest(folder/'moments.npz')!=p['geometry_sha256']:raise ValueError('Joint measurements changed')
    seed=p['seeds'][name]
    if digest(seed['path'])!=seed['sha256']:raise ValueError('Joint initialization changed')
    own=np.asarray(json.loads(Path(seed['path']).read_text())['model']['theta'])
    features,norm=model.response_features(a,train,p['feature_normalization']);_,field,nu=[x for x in p['configs'] if x[0]==name][0]
    dest=directory/name;dest.mkdir(exist_ok=True)
    if (dest/'freeze.json').exists():raise FileExistsError('Joint model already frozen')
    obj=model.JointADC(a,g,train,features,adc,penalty=p['penalty'],nodes=p['nodes'],time_field=field,rigidity_nu=nu,threads=4)
    initial=[model.parameters_initial(features,adc,train,1.,own),model.parameters_initial(features,adc,train,3.)];fits=[]
    for i,point in enumerate(initial):
        path=dest/('start_%d.json'%i)
        if path.exists():result=json.loads(path.read_text())
        else:
            rounds=[]
            for j in range(3):
                result=solve(obj,point,p['bounds']);rounds.append({k:v for k,v in result.items() if k!='theta'})
                print('JOINT_ADC_FIT',dataset,name,i,j,result['penalized_nll'],result['projected_gradient_max'],result['elapsed_seconds'],flush=True)
                if result['success'] and result['projected_gradient_max']<1e-3:break
                point=np.asarray(result['theta'])
            result['rounds']=rounds;write(path,result)
        fits.append(result)
    best=min(fits,key=lambda v:v['penalized_nll'])
    if not best['success'] or best['projected_gradient_max']>=1e-3:raise RuntimeError('Joint fit not converged')
    write(dest/'freeze.json',dict(complete=True,created_utc=stamp(),model=best,starts=fits,protocol_sha256=digest(directory/'protocol.json'),MC_truth_read=False,RICH_read=False))
    theta=np.asarray(best['theta']);ev=model.JointADC(a,g,test,features,adc,penalty=0,nodes=p['evaluation_nodes'],time_field=field,rigidity_nu=nu,threads=4)
    loss,grad,per=ev.evaluate(theta,events=True);out=model.reconstruct(theta,g['x'][test],ev.y,features[test])
    data={k:a[k][test] for k in ['run','event','R']}
    for k in ['joint_u','time_u','adc_u','weights','time_weights','adc_weights','time_injection_gain','adc_injection_gain','information','time_information','additional_adc_information']:data[k]=out[k]
    corrections={}
    for label,weight in [('time',out['time_weights']),('joint',out['weights'][:,:3].dot(base.D))]:
        solved=transport.solve_center(out[label+'_u'],a['R'][test],weight,g['moment'][test],theta[:2])
        if not solved['converged'].all():raise RuntimeError('Joint mean-loss inverse did not converge')
        data[label+'_transport_u']=solved['u'];corrections[label]=dict(max_inverse_residual=float(np.max(abs(solved['residual']))),median_correction=float(np.median(solved['correction'])),max_abs_correction=float(np.max(abs(solved['correction']))))
    b=np.linspace(-.1,.1,len(ev.y));changed=model.reconstruct(theta,g['x'][test],ev.y+b[:,None]*out['direction'],features[test]);injection=changed['joint_u']-out['joint_u']-b
    if np.max(abs(injection))>2e-12:raise ValueError('Joint coherent response failed')
    ev128=model.JointADC(a,g,test,features,adc,penalty=0,nodes=128,time_field=field,rigidity_nu=nu,threads=4);v128,_,p128=ev128.evaluate(theta,events=True)
    data['nll']=per;data['nll128']=p128;np.savez_compressed(dest/'predictions_before_truth.npz',**data)
    write(dest/'measurement_evaluation.json',dict(complete=True,created_utc=stamp(),events=len(per),mean_joint_nll=loss,
        quadrature_mean128_minus256=v128-loss,quadrature_rms128_minus256=float(np.sqrt(np.mean((p128-per)**2))),quadrature_max_abs128_minus256=float(np.max(abs(p128-per))),
        coherent_joint_injection_max_error=float(np.max(abs(injection))),time_only_injection_gain_percentiles=np.percentile(out['time_injection_gain'],[0,5,50,95,100]),
        model_information_width_gain_percentiles=np.percentile(1-np.sqrt(out['time_information']/out['information']),[0,5,50,95,100]),
        model_information_scope='Conditional linear-response, correctly specified common noise scale only; not an observed resolution gain',
        ADC_lambda=np.exp(theta[35:39]),material_kappa=theta[:2],loss_inverse=corrections,
        beta_above_one_retained={k:int((v<1).sum()) for k,v in data.items() if k.endswith('_u')},
        freeze_sha256=digest(dest/'freeze.json'),output_sha256=digest(dest/'predictions_before_truth.npz'),MC_truth_read=False,RICH_read=False))
    print('JOINT_ADC_FROZEN',dataset,name,'NLL',loss,'lambda',np.exp(theta[35:39]),'timegain',np.median(out['time_injection_gain']),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','fit']);parser.add_argument('--dataset',choices=['mc','flight'],required=True);parser.add_argument('--name',default='native_gaussianR');args=parser.parse_args()
    if args.stage=='prepare':prepare(args.dataset)
    else:fit(args.dataset,args.name)
