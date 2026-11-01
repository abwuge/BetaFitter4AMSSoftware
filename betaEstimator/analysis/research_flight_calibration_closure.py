#!/usr/bin/env python3
"""Measurement-only closure of the linearized calibration primitive.

Two synthetic detector domains have different calibration and drift. Each
fit sees measured arrays only. This is a numerical test, not AMS MC or a
claim about flight beta resolution.
"""
import json
import hashlib
from pathlib import Path
import numpy as np
import flight_calibration_core as core

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/flight_calibration_research'
PERIODS=6
P=4+3*PERIODS
PROTOCOL=dict(seed=21471,train_events_per_domain=50000,test_events_per_domain=12000,
              periods=PERIODS,loss_columns=4,clock_contrasts_per_period=3,
              sigmaA_ns=.045,sigmaB_ns=.065,
              scope='Exact linearized detector model; not AMS simulation',
              restriction='No MC-trained coefficients or truths passed to calibration',
              limitations='Prescribed correct template family and independent homoscedastic noise. This does not validate physical transport or selection modeling.')


def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parameters(domain):
    rng=np.random.RandomState(PROTOCOL['seed']+domain*101)
    loss=np.array([.035,.028,.032,.043])*(1. if domain==0 else 1.55)
    period=np.arange(PERIODS)
    offsets=np.column_stack([.03*np.sin(period+domain),.025*np.cos(period*.8+domain),
                             .01*(period-2.5)*(1. if domain==0 else -1.)])
    offsets+=rng.normal(0,.008,(PERIODS,3))
    return np.r_[loss,offsets.ravel()]


def generate(domain,n,seed):
    rng=np.random.RandomState(seed)
    x=np.array([-2.17,-2.04,2.04,2.17])[None,:]*rng.uniform(1.,1.3,(n,1))
    x+=rng.normal(0,.018,(n,4));x-=x.mean(axis=1)[:,None]
    epoch=rng.randint(0,PERIODS,n)
    high=rng.uniform(size=n)<.3
    sensitivity=np.where(high,.001,1.)
    amplitude=np.clip(rng.normal(0,1,(n,4)),-2.,2.)
    position=np.stack([.5*(x[:,0]+x[:,1]),.5*x[:,1],.5*x[:,2],.5*(x[:,2]+x[:,3])],axis=1)
    b=np.zeros((n,4,P))
    for k in range(4):
        b[:,:,k]=core.impulse_kernel(x,position[:,k])*(sensitivity*(1+.3*amplitude[:,k]))[:,None]
    for k in range(3):
        b[np.arange(n),k,4+3*epoch+k]=1.
        b[np.arange(n),3,4+3*epoch+k]=-1.
    u=np.where(high,rng.uniform(1.0002,1.003,n),rng.uniform(1.15,1.65,n))
    t0=rng.normal(0,1,n)
    sigma=PROTOCOL['sigmaA_ns'] if domain==0 else PROTOCOL['sigmaB_ns']
    theta=parameters(domain)
    times=t0[:,None]+u[:,None]*x+np.einsum('nmp,p->nm',b,theta)+rng.normal(0,sigma,(n,4))
    measured=dict(x=x,times=times,design=b,period=epoch,high=high)
    truth=dict(u=u,theta=theta,sigma=sigma)
    return measured,truth


def train(m):
    assert set(m)=={'x','times','design','period','high'}
    model=core.solve_global(m['x'],m['times'],m['design'],1.)
    sigma=np.sqrt(model['residual_sum_squares']/model['residual_dof'])
    model['estimated_sigma_ns']=float(sigma)
    return model


def prediction(m,theta):
    return core.reconstruct(m['x'],m['times'],m['design'],1.,theta)['u']


def stats(v):
    q=np.percentile(v,[16,84])
    return dict(mean=float(v.mean()),w68=float((q[1]-q[0])/2),rmse=float(np.sqrt(np.mean(v*v))))


def main():
    OUT.mkdir(exist_ok=True)
    (OUT/'closure_protocol.json').write_text(json.dumps(PROTOCOL,indent=2)+'\n')
    records=[];models=[];test_data=[];test_truth=[]
    for domain in [0,1]:
        m,tr=generate(domain,PROTOCOL['train_events_per_domain'],PROTOCOL['seed']+1000+domain)
        v,truth=generate(domain,PROTOCOL['test_events_per_domain'],PROTOCOL['seed']+2000+domain)
        model=train(m)
        models.append(model);test_data.append(v);test_truth.append(truth)
        # Freeze before using generator truths in any score.
        frozen=dict(theta=model['theta'].tolist(),estimated_sigma_ns=model['estimated_sigma_ns'],
                    data_rank=model['data_rank'],core_sha256=digest(core.__file__),protocol_sha256=digest(OUT/'closure_protocol.json'))
        (OUT/('closure_domain_%d_freeze.json'%domain)).write_text(json.dumps(frozen,indent=2)+'\n')
        pred=prediction(v,model['theta']);raw=prediction(v,np.zeros(P))
        oracle=prediction(v,truth['theta'])
        error=model['theta']-tr['theta']
        se=np.sqrt(np.diag(model['covariance_on_estimable_subspace']))*model['estimated_sigma_ns']
        eps=.007
        shift=dict(v,times=v['times']+eps*v['x']+2.31)
        response=np.max(np.abs(prediction(shift,model['theta'])-pred-eps))
        rawdesign=m['design'][:,:,:4]
        clock=np.zeros((len(m['x']),4,3))
        for k in range(3):clock[:,k,k]=1;clock[:,3,k]=-1
        pooled_design=np.concatenate([rawdesign,clock],axis=2)
        pooled=core.solve_global(m['x'],m['times'],pooled_design,1.)
        vclock=np.broadcast_to(clock[0],(len(v['x']),4,3))
        vp=np.concatenate([v['design'][:,:,:4],vclock],axis=2)
        pooled_u=core.reconstruct(v['x'],v['times'],vp,1.,pooled['theta'])['u']
        record=dict(domain=domain,train_n=len(m['x']),test_n=len(v['x']),rank=model['data_rank'],
                    estimated_sigma_ns=model['estimated_sigma_ns'],generator_sigma_ns=tr['sigma'],
                    max_abs_parameter_pull=float(np.max(np.abs(error/se))),
                    max_slope_and_t0_equivariance_error=float(response),
                    loss_parameters=model['theta'][:4].tolist(),loss_truth=tr['theta'][:4].tolist(),loss_standard_error=se[:4].tolist(),
                    metrics={k:stats(p-truth['u']) for k,p in [('raw',raw),('time_local',pred),('time_pooled',pooled_u),('known_generator_parameters',oracle)]},
                    high_speed_measured_u_below_one_fraction=float(np.mean(pred[v['high']]<1)))
        records.append(record)
    # Explicitly demonstrate why coefficients do not transfer across domains.
    transfer=prediction(test_data[1],models[0]['theta'])
    records[1]['domain_A_parameters_on_B']=stats(transfer-test_truth[1]['u'])
    # A global slope-shaped parameter is exactly invisible to self-consistency.
    n=120;x=np.tile(np.array([-2.17,-2.04,2.04,2.17]),(n,1))
    b=x[:,:,None];t=1.2*x
    caught=False
    try:core.solve_global(x,t,b,1.)
    except ValueError:caught=True
    hidden=core.solve_global(x,t,b,1.,allow_rank_deficient=True,penalty=np.eye(1))
    assert caught and hidden['data_rank']==0
    # Mean-zero shrinkage changes no ensemble mean, yet destroys unit response.
    rng=np.random.RandomState(910)
    u=1.2+rng.uniform(-.04,.04,1000)
    shrunk=u.mean()+.2*(u-u.mean())
    shrink=dict(mean_difference=float(shrunk.mean()-u.mean()),response=.2)
    # Two bounded latent decompositions give exactly the same measured times.
    # No energy sensor or mass-rigidity anchor is assumed in this example.
    n=2000;x=np.tile(np.array([-2.17,-2.04,2.04,2.17]),(n,1))
    z1=rng.uniform(-.01,.01,n);z2=rng.uniform(-.025,.025,n)
    ua=1.3+z1;ub=1.3+z2;noise=rng.normal(0,.05,x.shape)
    ta=(ua+z2)[:,None]*x+noise;tb=(ub+z1)[:,None]*x+noise
    twins=dict(max_measured_time_difference=float(np.max(np.abs(ta-tb))),
               true_u_difference_rmse=float(np.sqrt(np.mean((ua-ub)**2))),
               assumptions='Unknown slope-correlated noise allowed; no independent speed or known-noise anchor. Not an impossibility proof under stronger physical constraints.')
    result=dict(complete=True,protocol=PROTOCOL,domains=records,
                exact_slope_mode_rejected=caught,data_rank_with_slope_prior=hidden['data_rank'],
                mean_zero_shrink_counterexample=shrink,observational_twins=twins,
                core_sha256=digest(core.__file__),script_sha256=digest(__file__))
    (OUT/'linearized_closure.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
