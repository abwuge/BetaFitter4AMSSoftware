#!/usr/bin/env python3
"""Center-plane beta width, bias and gain from newly frozen ratio predictions."""
import json
import multiprocessing as mp
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import time_local_ratio as model
import evaluate_tofcalib_validation as binning

OUT = model.OUT
ROOT = model.ROOT
EV = model.local.w.ev


def summarize(task):
    rad,axis,low,high,values,old,new,reference,runs = task
    n = len(old)
    row = dict(low=float(low),high=float(high),n=n)
    if n >= 20:
        row.update(EV.compare(old,new,reference))
        row['x']=float(np.median(values))
        row['bootstrap']=EV.bootstrap(old,new,reference,runs,repeats=400)
    return rad,axis,row


def main(workers):
    freeze=json.loads((OUT/'freeze.json').read_text());verified=json.loads((OUT/'verification.json').read_text())
    if not freeze['complete'] or not verified['passed']:raise RuntimeError('Train/freeze/verify first')
    if verified['prediction_sha256']!=EV.digest(OUT/'predictions.npz'):raise RuntimeError('Prediction snapshot changed')
    p=model.local.load(OUT/'predictions.npz')
    with np.errstate(divide='raise',invalid='raise'):
        reference=1/(1/p['rich_beta'].astype(float)-p['delay_u_C12'])
        linear=1/p['raw_u'].astype(float)
        new=1/p['ratio_corrected_u']
    protocol=json.loads((OUT/'protocol.json').read_text())
    r_edges=np.unique(np.r_[binning.R_EDGES[binning.R_EDGES<protocol['coverage_GV'][1]],protocol['coverage_GV'][1]])
    tasks=[];overall={};old_neighbour_comparison={};groups={name:{axis:[] for axis in ['beta','rigidity']} for name in ['NaF','AGL']}
    for naf,rad in [(True,'NaF'),(False,'AGL')]:
        s=p['reference_mask']&(p['NaF'].astype(bool)==naf)
        if not np.isfinite(new[s]).all() or not np.isfinite(reference[s]).all():raise RuntimeError('Non-finite complete-cohort outputs')
        overall[rad]=EV.compare(linear[s],new[s],reference[s])
        old_neighbour_comparison[rad]=EV.compare(1/p['local_corrected_u'][s],new[s],reference[s])
        for axis,values,edges in [('beta',p['rich_beta'],binning.BETA_EDGES),('rigidity',p['R'],r_edges)]:
            for low,high in zip(edges[:-1],edges[1:]):
                mask=s&(values>=low)&(values<high)
                tasks.append((rad,axis,low,high,values[mask],linear[mask],new[mask],reference[mask],p['run'][mask]))
    with mp.Pool(workers) as pool:
        for rad,axis,row in pool.imap(summarize,tasks):
            groups[rad][axis].append(row)
    result=dict(complete=True,prediction_sha256=verified['prediction_sha256'],
        model_freeze_sha256=EV.digest(OUT/'freeze.json'),groups=groups,overall=overall,
        old_neighbour_same_event_comparison=old_neighbour_comparison,
        residual='1/u_TOF - 1/(1/beta_RICH - delay_u_C12)',reference_plane='L5',
        x_axis='Measured RICH beta (unchanged from original comparison) or measured R',
        uncertainty='400 paired whole-run bootstrap resamples, fixed trained models; no ratio-error/systematic propagation',
        evaluation_sample='Previously used run-separated flight validation events; newly generated predictions, not a new blind sample',
        no_output_clipping=True,NaF_AGL_separate=True)
    EV.dump(OUT/'center_results.json',result)
    dest=OUT/'figures';dest.mkdir(exist_ok=True)
    for rad in ['NaF','AGL']:
        folder=dest/rad;folder.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='ratio-trends-',dir='/tmp') as scratch:
            source=Path(scratch)/'points.txt'
            with source.open('w') as stream:
                stream.write('2\n')
                for axis in ['beta','rigidity']:
                    rows=[r for r in groups[rad][axis] if r['n']>=20]
                    bounds=(binning.BETA_EDGES[0],binning.BETA_EDGES[-1]) if axis=='beta' else tuple(protocol['coverage_GV'])
                    if rad=='AGL' and axis=='beta':bounds=(.955,1.02)
                    stream.write('%s %d %.17g %.17g\n'%(axis,len(rows),bounds[0],bounds[1]))
                    for row in rows:
                        ci=row['bootstrap']['beta_gain_pct_q025_q50_q975']
                        v=[row['x']]+[row[k]['w68']*1000 for k in ['baseline_beta','new_beta']]+[row[k]['mean']*1000 for k in ['baseline_beta','new_beta']]+[row['gain_beta_pct'],ci[0],ci[2]]
                        stream.write(' '.join('%.17g'%x for x in v)+'\n')
            subprocess.run(['root','-l','-b','-q','%s("%s","%s",true,"%s")'%(ROOT/'macro/plotTimeLocalTrends.C',source,folder,rad)],check=True)
    pages=[dest/rad/(axis+'_'+metric+'.pdf') for rad in ['NaF','AGL'] for axis in ['beta','rigidity'] for metric in ['width','bias','gain']]
    subprocess.run(['pdfunite']+[str(x) for x in pages]+[str(dest/'center_comparison.pdf')],check=True)
    print('CENTER_COMPARISON_COMPLETE',json.dumps(overall),flush=True)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=8)
    main(parser.parse_args().workers)
