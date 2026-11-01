#!/usr/bin/env python3
"""Recompute the MC check after ratio training, on the center beta plane."""
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import uproot
import time_local_ratio as model
import finish_material_comparison as material
import evaluate_tofcalib_validation as bins


def main(reuse_prepared=False):
    root=model.ROOT;out=model.OUT;ev=model.local.w.ev
    verified=json.loads((out/'verification.json').read_text())
    if not verified['passed']:raise RuntimeError('Verify new predictions before plotting')
    dest=out/'material_closure';dest.mkdir(exist_ok=True)
    material.OUT=dest
    if reuse_prepared:
        cached=json.loads((dest/'metrics.json').read_text())
        if cached['flight_sha256']!=verified['prediction_sha256']:
            raise RuntimeError('Prepared MC check is not associated with the new predictions')
    else:
        material.prepare(out/'predictions.npz')
    f=uproot.open(str(dest/'comparisons.root'))
    groups=[];results={}
    for rad in ['NaF','AGL']:
        a=f['mc_'+rad].arrays(library='np')
        residual=1/a['u_reference']-a['beta_mc']
        edges=bins.BETA_EDGES[bins.BETA_EDGES>=.955] if rad=='AGL' else bins.BETA_EDGES
        isotope_rows=[];isotopes={}
        for label in [0,1]:
            rows=[]
            for lo,hi in zip(edges[:-1],edges[1:]):
                s=(a['beta_mc']>=lo)&(a['beta_mc']<hi)&(a['isotope']==label);n=int(s.sum())
                if n<15:continue
                x=residual[s];source=a['source'][s]
                indices=[np.flatnonzero(source==v) for v in np.unique(source)]
                random=np.random.RandomState(20261016);boot=[]
                for _ in range(400):
                    ix=np.concatenate([indices[j] for j in random.randint(len(indices),size=len(indices))])
                    v=x[ix];boot.append([float(v.mean()),ev.width(v)])
                err=np.std(boot,axis=0,ddof=1)
                rows.append(dict(low=float(lo),high=float(hi),n=n,x=float(np.median(a['beta_mc'][s])),
                    mean=float(x.mean()),w68=ev.width(x),mean_error=float(err[0]),w68_error=float(err[1])))
            isotope_rows.append(rows)
            isotopes['C%d'%(12+label)]=dict(n=int((a['isotope']==label).sum()),rows=rows,
                overall_beta=ev.stats(residual[a['isotope']==label]))
        results[rad]=dict(n=len(residual),overall_beta=ev.stats(residual),isotopes=isotopes)
        for metric in ['mean','w68']:
            groups.append((rad,metric,isotope_rows))
    ev.dump(dest/'center_closure.json',dict(complete=True,results=results,
        residual='1/(1/beta_RICH - delay_u_C12) - beta_MC(L5)',
        material_model='Same fixed nominal C12 propagation as flight validation',
        errors='Standard deviation of 400 whole-source-file bootstrap replicates',
        sample='Existing 2.5--6 GV MC reference cohort, material evaluation freshly recomputed after new-model freeze',
        model_dependence='Material/RICH closure is independent of the trained TOF correction'))
    with tempfile.TemporaryDirectory(prefix='ratio-closure-',dir='/tmp') as scratch:
        source=Path(scratch)/'points.txt'
        with source.open('w') as stream:
            stream.write(str(len(groups))+'\n')
            for rad,metric,isotope_rows in groups:
                stream.write('%s %s %.17g %.17g %d %.17g\n'%(rad,metric,.78 if rad=='NaF' else .955,1.02,results[rad]['n'],results[rad]['overall_beta']['mean']*1000))
                for rows in isotope_rows:
                    stream.write(str(len(rows))+'\n')
                    for row in rows:stream.write('%.17g %.17g %.17g %d\n'%(row['x'],row[metric]*1000,row[metric+'_error']*1000,row['n']))
        subprocess.run(['root','-l','-b','-q','%s("%s","%s")'%(root/'macro/plotRatioMaterialClosure.C',source,dest)],check=True)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--reuse-prepared',action='store_true')
    main(parser.parse_args().reuse_prepared)
