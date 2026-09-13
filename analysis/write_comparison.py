#!/usr/bin/env python3
"""Write aligned Hybrid/new-estimator comparison trees using system uproot."""
import argparse, json, os
import numpy as np
import uproot

CHARGES=(1,2,6,8)

def load(path):
    with np.load(path,allow_pickle=False) as f:return {k:f[k] for k in f.files}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cache',default='results/comparison/cache');ap.add_argument('--input',default='results/comparison');ap.add_argument('--output',default='results/comparison');args=ap.parse_args()
    if not os.path.isdir(args.output):os.makedirs(args.output)
    summary={'charges':{}}
    for z in CHARGES:
        d=load(os.path.join(args.cache,'validation_Z%d.npz'%z));u=np.load(os.path.join(args.input,'new_inverse_beta_Z%d.npy'%z));new_beta=1/u
        common=(d['charge']==z)&(d['mch'].astype(np.int32)==z)&np.isfinite(d['mcBeta'])&(d['mcBeta']>=.25)&(d['mcBeta']<=1)
        for k in ['linearBeta','nonlinearBeta','hybridBeta']:common &= np.isfinite(d[k])&(d[k]>0)
        common &= np.isfinite(new_beta)&(new_beta>0)
        dest=os.path.join(args.output,'comparison_Z%d.root'%z)
        with uproot.recreate(dest) as f:
            f['comparison']={
                'run':d['run'][common].astype(np.uint32),'event':d['event'][common].astype(np.uint32),'sourceId':d['source_id'][common].astype(np.int32),'Z':d['charge'][common].astype(np.int32),
                'mcBeta':d['mcBeta'][common].astype(np.float64),'hybridBeta':d['hybridBeta'][common].astype(np.float64),'newBeta':new_beta[common].astype(np.float64),
                'linearBeta':d['linearBeta'][common].astype(np.float64),'nonlinearBeta':d['nonlinearBeta'][common].astype(np.float64),
            }
        n06=int(np.count_nonzero(common&(d['mcBeta']>=.6)))
        summary['charges'][str(z)]={'common_entries':int(common.sum()),'beta_ge_06_entries':n06,'output':dest}
        print('Z=%d common=%d beta>=0.6=%d'%(z,common.sum(),n06),flush=True)
    with open(os.path.join(args.output,'comparison_summary.json'),'w') as f:json.dump(summary,f,indent=2)
if __name__=='__main__':main()
