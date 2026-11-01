"""Free timing weights with a continuous-mass, first-order L5 transport inverse.

The event mass is p*sqrt(u*u-1), with measured p=Z*R. No isotope probability
or source-spectrum value is evaluated. The forward extension is zero loss
for u<=1 so that measurement outputs above beta=1 remain representable.
"""
import numpy as np
import flight_mean_loss as physical
from research_conditional_output import weights


def compress(moment):
    m=np.asarray(moment,float)
    if m.ndim!=4 or m.shape[1:]!=(4,2,9) or not np.isfinite(m).all():
        raise ValueError('Expected finite n x 4 x 2 x 9 material moments')
    return np.sum(m*physical.ZA,axis=-1),np.sum(m*physical.ZA*np.log(physical.I_MEV),axis=-1)


def continuous_kernel(u,rigidity,moment0,moment1,charge=6):
    """First-order times by layer and material group, in ns."""
    u=np.asarray(u,float);r=np.asarray(rigidity,float)
    if u.ndim!=1 or r.shape!=u.shape or not np.isfinite(u).all() or not np.isfinite(r).all() or np.any(r<=0):
        raise ValueError('Finite inverse speeds and positive measured rigidities required')
    mm0=np.asarray(moment0,float);mm1=np.asarray(moment1,float)
    if mm0.shape!=(len(u),4,2) or mm1.shape!=mm0.shape:raise ValueError('Compressed moment dimensions differ')
    out=np.zeros_like(mm0);valid=u>1
    if not valid.any():return out
    v=u[valid];p=charge*r[valid];d=(v-1)*(v+1);b2=1/(v*v);g2=v*v/d
    mass=p*np.sqrt(d);me=.511;mass_mev=1000*mass
    tmax=2*me*b2*g2/(1+2*np.sqrt(g2)*me/mass_mev+(me/mass_mev)**2)
    coefficient=d/p*(.000307075*charge*charge/b2)
    bracket=.5*np.log(2*me*b2*g2*tmax)-b2
    out[valid]=coefficient[:,None,None]*(mm0[valid]*bracket[:,None,None]-mm1[valid])
    return out


def solve_center(slope,rigidity,weight,moment,kappa,tolerance=2e-12,maxiter=40):
    """Invert slope = u_L5 + w dot B(u_L5, measured R) kappa.

    Newton steps are capped for numerical stability, not to constrain beta.
    The returned convergence mask must be checked by callers.
    """
    slope=np.asarray(slope,float);r=np.asarray(rigidity,float);w=np.asarray(weight,float);k=np.asarray(kappa,float)
    if slope.ndim!=1 or w.shape!=(len(slope),4) or r.shape!=slope.shape:raise ValueError('Incompatible timing inputs')
    if k.shape!=(2,) or not np.isfinite(k).all() or np.any(k<0):raise ValueError('Two nonnegative material scales required')
    if not np.isfinite(slope).all() or not np.isfinite(w).all():raise ValueError('Nonfinite time slope or weights')
    m0,m1=compress(moment)
    def delta(v):return np.einsum('nl,nlg,g->n',w,continuous_kernel(v,r,m0,m1),k)
    u=slope.copy();iterations=0
    for i in range(maxiter):
        value=u+delta(u)-slope;iterations=i+1
        if np.max(np.abs(value))<tolerance:break
        step=1e-5*np.maximum(1,np.abs(u))
        derivative=1+(delta(u+step)-delta(u-step))/(2*step)
        if np.any(~np.isfinite(derivative)) or np.any(derivative<=.25):
            raise FloatingPointError('Transport inverse is not locally monotone in this state')
        u-=np.clip(value/derivative,-.1,.1)
    residual=u+delta(u)-slope
    return dict(u=u,correction=slope-u,residual=residual,converged=np.abs(residual)<tolerance,iterations=iterations,
                forward_extension='Zero material correction at u<=1; no clipping of measured inverse speed')


def center_from_times(theta,x,t,rigidity,moment,kappa=None):
    w=weights(np.asarray(theta,float),np.asarray(x,float));s=np.sum(w*np.asarray(t,float),axis=1)
    scale=np.asarray(theta[:2] if kappa is None else kappa,float)
    return solve_center(s,rigidity,w,moment,scale)
