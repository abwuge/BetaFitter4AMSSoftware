"""Necessary existing RICH cuts, before extracting missing reference fields.

This is evaluation-only. Missing track matching and beta-error cuts remain
mandatory in the final frozen reference selection; no residual is used here.
"""
import numpy as np


def reference_preselection(a):
    b=np.asarray(a['rich_beta'][:,0],float)
    naf=a['rich_NaF'].astype(bool)
    xy=np.max(np.abs(a['rich_pos'][:,:2]),axis=1)
    rad=np.sum(a['rich_pos'][:,:2]**2,axis=1)
    npe=a['rich_npe'];bc=a['rich_BetaConsistency'];pc=a['rich_PMTChargeConsistency']
    ok=(a['irich']>=0)&(a['itrtrack']>=0)
    ok &= a['rich_good'].astype(bool)&a['rich_clean'].astype(bool)
    ok &= (a['rich_pb']>=.01)&(a['rich_pmt']>=3)&(a['rich_used']>=4)
    ok &= np.isfinite(b)&(b>=.78)&(b<=1.02)&np.where(naf,b>=.78,b>=.955)
    ok &= np.isfinite(bc)&(bc>=0)&(bc<=np.where(naf,.01,.005))
    ok &= np.isfinite(pc)&(pc>=0)&(pc<=15)
    ok &= np.all(np.isfinite(npe),axis=1)&(npe[:,2]>0)
    ok &= (npe[:,0]>=.4*npe[:,2])&(npe[:,1]>=np.where(naf,1.,2.))
    ok &= np.where(naf,xy<=17,(xy>=19)&(rad<=3500)&~np.isin(a['rich_tile'],[3,7,12,20,87,100,108]))
    ok &= np.all(np.isfinite(a['rich_pos']),axis=1)
    return ok & (a['rich_q'][:,1]>4)&(a['rich_q'][:,1]<8)
