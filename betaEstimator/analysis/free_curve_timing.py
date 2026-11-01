"""Refine a continuous time/ADC output against its material-aware mean curve.

This minimizes a fixed-model Mahalanobis residual. It does not evaluate an
isotope or source posterior and does not claim a global nonlinear optimum.
"""
import numpy as np
import flight_conditional_likelihood as timing
import free_timing_transport as material


def mean_curve(u, rigidity, direction, moment0, moment1, kappa):
    loss=np.einsum('nlg,g->nl',material.continuous_kernel(u,rigidity,moment0,moment1),kappa)
    out=direction*np.asarray(u)[:,None]
    out[:,:3]+=loss.dot(timing.D.T)
    return out


def refine(initial, centered, direction, precision, rigidity, moment, kappa, tolerance=2e-10, max_iterations=30):
    initial=np.asarray(initial,float);y=np.asarray(centered,float);d=np.asarray(direction,float)
    q=np.asarray(precision,float);r=np.asarray(rigidity,float);k=np.asarray(kappa,float)
    if initial.ndim!=1 or y.shape!=(len(initial),7) or d.shape!=y.shape or q.shape!=(7,7) or r.shape!=initial.shape or k.shape!=(2,):
        raise ValueError('Incompatible continuous curve inputs')
    if not all(np.isfinite(value).all() for value in [initial,y,d,q,r,k]) or np.any(k<0) or np.linalg.eigvalsh(q)[0]<=0:
        raise ValueError('Nonfinite or nonphysical curve inputs')
    m0,m1=material.compress(moment)
    def curve(u):return mean_curve(u,r,d,m0,m1,k)
    def cost(u):
        residual=y-curve(u)
        return np.einsum('ni,ij,nj->n',residual,q,residual)
    def derivative(u):
        step=1e-5*np.maximum(1,abs(u))
        return (curve(u+step)-curve(u-step))/(2*step[:,None])
    u=initial.copy();before=cost(u);current=before.copy();accepted_steps=np.zeros(len(u),int)
    converged=np.zeros(len(u),bool);line_failures=np.zeros(len(u),bool)
    for iteration in range(max_iterations):
        residual=y-curve(u);jacobian=derivative(u)
        information=np.einsum('ni,ij,nj->n',jacobian,q,jacobian)
        score=np.einsum('ni,ij,nj->n',jacobian,q,residual)
        step=score/information
        converged=abs(step)<=tolerance
        active=~converged
        if not active.any():break
        step=np.clip(step,-.1,.1)
        scale=np.ones(len(u));chosen=np.zeros(len(u),bool);proposal=u.copy();next_cost=current.copy()
        for line in range(25):
            waiting=active&~chosen
            if not waiting.any():break
            trial=u+scale*step;value=cost(trial)
            take=waiting&(value <= current - 1e-4*2*scale*score*step + 1e-13*np.maximum(1,current))
            proposal[take]=trial[take];next_cost[take]=value[take];chosen|=take
            scale[waiting&~take]*=.5
        line_failures|=active&~chosen
        accepted_steps+=chosen.astype(int);u=proposal;current=next_cost
        if not chosen.any():break
    residual=y-curve(u);jacobian=derivative(u)
    info=np.einsum('ni,ij,nj->n',jacobian,q,jacobian)
    stationary=np.einsum('ni,ij,nj->n',jacobian,q,residual)/info
    converged=abs(stationary)<=tolerance
    if not np.isfinite(u).all() or np.max(current-before)>1e-9*np.maximum(1,np.max(before)):
        raise FloatingPointError('Curve refinement increased the objective')
    return dict(u=u,converged=converged,score_in_u_units=stationary,initial_objective=before,objective=current,
        accepted_steps=accepted_steps,line_search_failure=line_failures,iterations=iteration+1,
        inverse_speed_change=u-initial)
