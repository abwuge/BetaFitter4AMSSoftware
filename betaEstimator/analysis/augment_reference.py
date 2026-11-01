#!/usr/bin/env python3
"""Add measured RICH quality fields to the sampled reference cache."""
import sys,json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np,uproot
FIELDS=['rich_pmt','rich_npe','rich_PMTChargeConsistency','rich_theta','rich_n','rich_cstat']
def read_one(rec):
    sid,path=rec
    with uproot.open(path) as f: a=f['amstreea'].arrays(FIELDS,library='np')
    return sid,a
def main():
    dest=Path(sys.argv[1])
    audit=json.loads((dest/'reference_audit.json').read_text())
    data={}
    index={}
    for fn in sorted(dest.glob('*_Z*.npz')):
        with np.load(str(fn),allow_pickle=False) as f:d={k:f[k] for k in f.files}
        for sid in np.unique(d['source_id']):
            index.setdefault(int(sid),[]).append((fn.name,np.flatnonzero(d['source_id']==sid)))
        data[fn.name]=d
    rec=[(s['source_id'],s['path']) for s in audit['sources']]
    with ProcessPoolExecutor(max_workers=3) as pool:
        for sid,a in pool.map(read_one,rec):
            for name,idx in index.get(sid,[]):
                d=data[name];entries=d['source_entry'][idx]
                for k,v in a.items():
                    if k not in d:d[k]=np.zeros((len(d['run']),)+v.shape[1:],dtype=v.dtype)
                    d[k][idx]=v[entries]
    for name,d in data.items():np.savez_compressed(str(dest/name),**d)
    audit['supplemental_measured_branches']=FIELDS
    (dest/'reference_audit.json').write_text(json.dumps(audit,indent=2))
    print('Supplemented measured RICH quality fields',flush=True)
if __name__=='__main__':main()
