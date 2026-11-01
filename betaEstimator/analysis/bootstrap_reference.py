#!/usr/bin/env python3
"""Paired source-block uncertainty for frozen predictions, never model selection."""
import argparse,json
from pathlib import Path
import numpy as np,uproot
def width(v):
    q=np.percentile(v,[16,84]);return float((q[1]-q[0])/2)
def weighted_width(v,group,counts):
    w=counts[group];cum=np.cumsum(w);total=cum[-1]
    qs=[]
    for p in (.16,.84):
        rank=p*(total-1);lo=int(np.floor(rank));hi=int(np.ceil(rank))
        a=v[np.searchsorted(cum,lo+1)];b=v[np.searchsorted(cum,hi+1)]
        qs.append(a+(b-a)*(rank-lo))
    return qs[1]-qs[0]
def compare(base,new,groups,seed,repeats=200):
    unique,gid=np.unique(groups,return_inverse=True);n=len(unique)
    a=np.argsort(base);b=np.argsort(new)
    av,ag=base[a],gid[a];bv,bg=new[b],gid[b]
    rng=np.random.RandomState(seed);samples=[]
    for i in range(repeats):
        counts=np.bincount(rng.randint(0,n,n),minlength=n)
        wa=weighted_width(av,ag,counts);wb=weighted_width(bv,bg,counts)
        samples.append(100*(1-wb/wa))
    return {'N':len(base),'sources':n,'baseline_width68':width(base),'new_width68':width(new),
        'narrowing_percent':100*(1-width(new)/width(base)),
        'paired_source_bootstrap_95':np.percentile(samples,[2.5,97.5]).tolist(),'replicates':repeats}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--sample',choices=['flight','historical','fresh'],required=True);args=ap.parse_args()
    out={'sample':args.sample,'uncertainty':'fixed-model source-file bootstrap; excludes training/model/systematic uncertainty','charges':{}}
    for z in (1,2,6,8):
        if args.sample=='flight':
            p=Path('results/reference_training')/('flight_Z%d.root'%z)
            with uproot.open(str(p)) as f:d=f['reference'].arrays(library='np')
            beta=d['richBeta'];base=d['linearU'];new=d['selectedU'];old=d['oldNllU'];ranges=[(.985,1.02),(.99,1.02)]
        else:
            p=Path('results/reference_comparison')/args.sample/('comparison_Z%d.root'%z)
            with uproot.open(str(p)) as f:d=f['comparison'].arrays(library='np')
            beta=d['mcBeta'];base=1/d['hybridBeta'];new=1/d['newBeta'];old=d['old_nllU'];ranges=[(.6,1.),(.95,1.)]
        out['charges'][str(z)]=[]
        for low,high in ranges:
            m=(beta>=low)&(beta<=high);r=1/beta[m]
            row={'beta_low':low,'beta_high':high}
            row['selected_vs_baseline']=compare(base[m]-r,new[m]-r,d['sourceId'][m],2181+z)
            row['selected_vs_old_nll']=compare(old[m]-r,new[m]-r,d['sourceId'][m],3181+z)
            out['charges'][str(z)].append(row)
    dest=Path('results/reference_training') if args.sample=='flight' else Path('results/reference_comparison')/args.sample
    (dest/'bootstrap_summary.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
