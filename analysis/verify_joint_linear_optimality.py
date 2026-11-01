#!/usr/bin/env python3
"""Check the exact limited optimality certificate of frozen linear estimators."""
import json
from pathlib import Path
import numpy as np
import joint_adc_likelihood as joint
import constrained_joint_adc as constrained
import periodic_joint_adc as period
from research_amplitude_noise_20261013 import inputs,CONFIGS
from research_flight_conditional import write,digest
from research_conditional_mc import stamp

DEST=period.OUT/'linear_certificate'


def certificate(q,direction,constraints,weights,rng):
    cov=np.linalg.inv(q);n=len(direction);response=np.zeros((n,constraints.shape[2]));response[:,0]=1
    null=rng.normal(size=(n,7))
    gram=np.einsum('nik,nil->nkl',constraints,constraints)
    coefficient=np.linalg.solve(gram,np.einsum('nik,ni->nk',constraints,null))
    null-=np.einsum('nik,nk->ni',constraints,coefficient)
    null*=.03
    original=np.einsum('ni,ij,nj->n',weights,cov,weights)
    moved=np.einsum('ni,ij,nj->n',weights+null,cov,weights+null)
    predicted=np.einsum('ni,ij,nj->n',null,cov,null)
    actual_response=np.einsum('nik,ni->nk',constraints,weights)
    null_error=float(np.max(abs(np.einsum('nik,ni->nk',constraints,null))))
    constraint_error=float(np.max(abs(actual_response-response)))
    gap_error=float(np.max(abs(moved-original-predicted)))
    if constraint_error>2e-10 or null_error>2e-10 or gap_error>1e-11 or np.any(moved-original < -1e-12):
        raise ValueError('Linear optimality identity failed')
    return dict(events=n,constraint_max_error=constraint_error,null_max_error=null_error,
        variance_gap_identity_max_error=gap_error,scale_variance_percentiles=np.percentile(original,[0,5,50,95,100]))


def main():
    DEST.mkdir(exist_ok=True)
    if (DEST/'verification.json').exists():raise FileExistsError('Linear certificate already saved')
    rng=np.random.RandomState(7126);rows=[]
    for dataset in ['flight','mc']:
        folder,a,g,adc,domain,train,test=inputs(dataset)
        selected=np.flatnonzero(test)[:128]
        original_norm=json.loads((joint.OUT/dataset/'protocol.json').read_text())['feature_normalization']
        features,_=joint.response_features(a,train,original_norm)
        x=g['x'][selected]
        for family,field,_ in CONFIGS[dataset]:
            source=joint.OUT/dataset/family/'freeze.json';theta=np.asarray(json.loads(source.read_text())['model']['theta'])
            observation=joint.observation({field:a[field][selected]},adc[selected],field)
            rec=joint.reconstruct(theta,x,observation,features[selected]);q=joint.precision(theta)[1];d=rec['direction']
            checks={'joint':certificate(q,d,d[:,:,None],rec['weights'],rng)}
            for variant in ['speed_null','speed_gain_null']:
                con=constrained.reconstruct(theta,x,observation,features[selected],variant)
                checks[variant]=certificate(q,d,con['constraint_matrix'],con['weights'],rng)
            variance_shift=[]
            for j in range(len(d)):
                shifted=rec['covariance']+.1*np.outer(d[j],d[j]);precision=np.linalg.inv(shifted)
                w=precision.dot(d[j]);w/=d[j].dot(w)
                variance_shift.append(np.max(abs(w-rec['weights'][j])))
            if max(variance_shift)>2e-10:raise ValueError('Signal-direction covariance invariance failed')
            rows.append(dict(dataset=dataset,family=family,model_sha256=digest(source),checks=checks,
                signal_direction_covariance_invariance_max_error=max(variance_shift),
                ADC_information_identity_max_error=float(np.max(abs(rec['information']-rec['time_information']-rec['additional_adc_information'])))))
    write(DEST/'verification.json',dict(complete=True,created_utc=stamp(),models=rows,source_sha256=digest(Path(__file__)),
        theorem='For fixed positive-definite covariance and fixed linear response constraints C^T w=b, w*=Q C (C^T Q C)^-1 b uniquely minimizes w^T Sigma w. Every feasible perturbation v has exact variance gap v^T Sigma v.',
        elliptical_extension='With a fixed linear mean d*u and centered multivariate Student noise with nu>2, the unconstrained GLS is also minimum squared-risk among translation-equivariant estimators: its error has conditional mean zero given the orthogonal residual.',
        boundaries=['This is a statement about a specified linear conditional response, not a global detector optimum.',
            'Estimated mean, covariance, era dependence and material parameters are not known truth.',
            'The physical material mean is nonlinear in u; applying a scalar material inverse does not make the seven-dimensional linear theorem a nonlinear global-optimality certificate.',
            'Feature expansion, different response laws, non-elliptical noise and calibration uncertainty remain outside the certified class.'],
        beta_references_read=False,production_acceptance=False))
    print('JOINT_LINEAR_CERTIFICATE',len(rows),'models verified; scope is fixed linear response only',flush=True)


if __name__=='__main__':main()
