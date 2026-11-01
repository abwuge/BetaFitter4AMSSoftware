"""Measured-likelihood curvature scaling, without changing model or penalty."""
import time
import numpy as np
from scipy.optimize import minimize


def projected_gradient(theta,gradient,bounds):
    p=np.asarray(gradient,float).copy()
    for i,(low,high) in enumerate(bounds):
        if theta[i]<=low+1e-7 and p[i]>0:p[i]=0
        if theta[i]>=high-1e-7 and p[i]<0:p[i]=0
    return p


def solve_scaled(objective,initial,bounds,maximum_iterations=6000):
    start=time.monotonic();theta=np.asarray(initial,float);bb=np.asarray(bounds,float)
    if bb.shape!=(len(theta),2) or np.any(theta<bb[:,0]-1e-8) or np.any(theta>bb[:,1]+1e-8):raise ValueError('Scaled optimizer starting bounds')
    diagonal=np.zeros(len(theta))
    for j in range(len(theta)):
        step=1e-4*max(1.,abs(theta[j]));up=theta.copy();down=theta.copy();up[j]+=step;down[j]-=step
        diagonal[j]=(objective.evaluate(up)[1][j]-objective.evaluate(down)[1][j])/(2*step)
    # Curvature is only a coordinate preconditioner. No Hessian uncertainty is inferred.
    scale=np.clip(1/np.sqrt(np.maximum(abs(diagonal),1e-6)),.001,10.)
    scaled_bounds=[((lo-v)/s,(hi-v)/s) for (lo,hi),v,s in zip(bounds,theta,scale)]
    trace=[];calls=[0]
    def fun(z):
        value,gradient=objective.evaluate(theta+scale*z);calls[0]+=1
        return value,gradient*scale
    def callback(z):
        if len(trace)%50==0:
            value,gradient=objective.evaluate(theta+scale*z)
            trace.append(dict(iteration=len(trace),value=float(value),projected_gradient_max=float(np.max(abs(projected_gradient(theta+scale*z,gradient,bounds))))))
        else:trace.append(None)
    result=minimize(fun,np.zeros(len(theta)),jac=True,method='L-BFGS-B',bounds=scaled_bounds,callback=callback,
        options=dict(maxiter=maximum_iterations,maxls=50,ftol=2e-13,gtol=1e-8,maxcor=80))
    fitted=theta+scale*result.x;value,gradient=objective.evaluate(fitted);unpenalized=objective.evaluate(fitted,regularized=False)[0]
    pg=projected_gradient(fitted,gradient,bounds)
    return dict(theta=fitted.tolist(),penalized_nll=float(value),nll=float(unpenalized),success=bool(result.success),message=str(result.message),
        iterations=int(result.nit),function_evaluations=int(result.nfev),gradient_max=float(np.max(abs(gradient))),projected_gradient_max=float(np.max(abs(pg))),
        elapsed_seconds=time.monotonic()-start,preconditioner_scale=scale.tolist(),observed_diagonal=diagonal.tolist(),
        trace=[row for row in trace if row is not None],optimizer='Curvature-scaled L-BFGS-B with memory80; same original model, bounds and penalty')
