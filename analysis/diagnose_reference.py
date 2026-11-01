#!/usr/bin/env python3
"""Variance decomposition of already-frozen calibration on final references."""
import json
from pathlib import Path
import numpy as np,uproot
def main():
    result={}
    for sample in ('flight','historical','fresh'):
        out={}
        for z in (1,2,6,8):
            path=Path('results/reference_training')/('flight_Z%d.root'%z) if sample=='flight' else Path('results/reference_comparison')/sample/('comparison_Z%d.root'%z)
            if not path.exists():continue
            with uproot.open(str(path)) as f:d=f['reference' if sample=='flight' else 'comparison'].arrays(library='np')
            b=d['richBeta'] if sample=='flight' else d['mcBeta']
            m=(b>=.99)&(b<=(1.02 if sample=='flight' else 1.))
            u0=d['linearU'][m] if sample=='flight' else d['olsU'][m]
            u=d['selectedU'][m]
            e=u0-1/b[m];c=u0-u
            if len(e)<20:continue
            ec=e-e.mean();cc=c-c.mean()
            ve=float(np.mean(ec*ec));vc=float(np.mean(cc*cc));cov=float(np.mean(ec*cc))
            after=float(np.var(e-c))
            out[str(z)]={'N':len(e),'variance_before':ve,'correction_variance':vc,
                'twice_error_correction_covariance':2*cov,'variance_after':after,
                'identity_error':after-(ve+vc-2*cov),'mean_before':float(e.mean()),
                'mean_correction':float(c.mean()),'mean_after':float((e-c).mean())}
        result[sample]=out
    Path('results/reference_comparison/variance_diagnostic.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
