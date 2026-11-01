"""Mean material loss and bidirectional flight-time integration at L5.

This module has no fitted coefficients, isotope abundance or detector labels.
It reproduces the elemental Bethe expression used by betaFitter. Density,
shell and nuclear-interaction corrections are not added by this module.
"""
import math
import numpy as np

C=29.9792458
ZA=np.array([.99212,.49954,.49976,.50002,.47372,.47847,.48181,.49848,.39575])
I_MEV=np.array([19.2e-6,78e-6,82e-6,95e-6,115e-6,149e-6,166e-6,173e-6,823e-6])
ATOMIC_MASS=np.array([1.0079,12.0107,14.0067,15.9994,18.9984,22.9897,26.9815,28.0855,207.2])

def stopping_power(beta,mass,charge=6):
    if not (np.isfinite(beta) and np.isfinite(mass) and 0<beta<1 and mass>0):raise ValueError('Physical beta and positive mass required')
    b2=beta*beta;g2=1/(1-b2);me=.511;m=mass*1000
    tmax=2*me*b2*g2/(1+2*math.sqrt(g2)*me/m+(me/m)**2)
    return .000307075*charge*charge*ZA/b2*(.5*np.log(2*me*b2*g2*tmax/I_MEV**2)-b2)

def inverse_beta(energy,mass):
    if not energy>mass:return np.nan
    return energy/math.sqrt((energy-mass)*(energy+mass))

def energy_loss(energy,mass,column,charge):
    u=inverse_beta(energy,mass)
    if not np.isfinite(u):return np.nan
    sp=stopping_power(1/u,mass,charge)
    return np.dot(sp,column)

def propagate_cells(u,mass,kappa,high,low,length,column,group,at_high,at_low,center_z,tof_z,substeps=1,charge=6):
    """Four times in ns; column is g/cm2 per element; loss is positive GeV."""
    if not (np.isfinite(u) and np.isfinite(mass) and u>1 and mass>0 and len(kappa)==2 and substeps>=1):raise ValueError('Invalid physical propagation input')
    if not np.isfinite(kappa).all() or np.any(np.asarray(kappa)<0):raise ValueError('Nonnegative finite loss scales required')
    energy0=mass*u/math.sqrt(u*u-1)
    times=np.full(4,np.nan);losses=np.zeros(2)
    n=len(high)
    for side in range(2):
        sign=1. if side==0 else -1.
        energy=energy0;time=0.
        for j in range(n):
            k=n-1-j if side==0 else j
            if side==0:
                if low[k]<center_z-1e-5 or high[k]>tof_z[0]+1e-5:continue
            else:
                if high[k]>center_z+1e-5 or low[k]<tof_z[3]-1e-5:continue
            col=column[k]*(kappa[group[k]]/substeps)
            for sub in range(substeps):
                first=energy_loss(energy,mass,col,charge)
                middle=energy+sign*.5*first
                loss=energy_loss(middle,mass,col,charge)
                next_energy=energy+sign*loss
                if not np.isfinite(next_energy) or next_energy<=mass:return np.full(4,np.nan),np.full(2,np.nan)
                time-=sign*.5*length[k]/(C*substeps)*(inverse_beta(energy,mass)+inverse_beta(next_energy,mass))
                losses[group[k]]+=loss
                energy=next_energy
            layer=at_high[k] if side==0 else at_low[k]
            if layer>=0:times[layer]=time
    return times,losses

def first_order(u,mass,moment,charge=6):
    if not (np.isfinite(u) and np.isfinite(mass) and u>1 and mass>0):raise ValueError('Physical inverse beta and positive mass required')
    p=mass/math.sqrt(u*u-1)
    sp=stopping_power(1/u,mass,charge)
    out=np.empty((4,2))
    for l in range(4):
        for g in range(2):out[l,g]=mass*mass/(p*p*p)*np.dot(moment[l,g],sp)
    return out
