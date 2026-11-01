#!/usr/bin/env python3
"""Prepare measured flight inputs and native geometry for a small mean-loss model."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import uproot

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/flight_calibration_research'
OUT=BASE/'mean_loss_likelihood'
C=29.9792458
ATOMIC_MASS=np.array([1.0079,12.0107,14.0067,15.9994,18.9984,22.9897,26.9815,28.0855,207.2])

def load(path):
    with np.load(path,allow_pickle=False) as f:return {k:f[k] for k in f.files}

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def keys(a):return (a['run'].astype(np.uint64)<<np.uint64(32))|a['event'].astype(np.uint64)

def match(a,b):
    x,y=keys(a),keys(b);o=np.argsort(y);j=np.searchsorted(y[o],x);safe=np.minimum(j,len(o)-1)
    return o[safe],(j<len(o))&(y[o[safe]]==x)

def prepare(limit):
    OUT.mkdir(exist_ok=True)
    a=load(BASE/'four_layer_timing_newcalib/measured.npz')
    with np.load(ROOT/'results/time_local_beta/training.npz') as f:b={k:f[k] for k in ['run','event','R','relR']}
    j,ok=match(a,b);ok&=(b['R'][j]>=2.5)&(b['R'][j]<6)
    ix=np.flatnonzero(ok)
    # Deterministic event sampling is independent of timing residuals and references.
    ranks=np.array([int.from_bytes(hashlib.sha256(('%d:%d'%(a['run'][i],a['event'][i])).encode()).digest()[:8],'big') for i in ix],np.uint64)
    ix=np.sort(ix[np.argsort(ranks)[:min(limit,len(ix))]])
    data={k:v[ix] for k,v in a.items()};data.update(R=b['R'][j[ix]],relR=b['relR'][j[ix]])
    n=load(BASE/'four_layer_timing_native/measured.npz');k,found=match(data,n)
    if not found.all():raise ValueError('Native identity missing')
    data['t_native']=n['t'][k]
    records=json.loads((ROOT/'results/tofcalib_validation/preparation.json').read_text())['records'];lookup={int(r['source']):r['path'] for r in records}
    data['tk_pos']=np.zeros((len(ix),9,3));data['tk_dir']=np.zeros((len(ix),9,3));source_notes=[]
    for num,source in enumerate(np.unique(data['source'])):
        pick=np.flatnonzero(data['source']==source)
        with uproot.open(lookup[int(source)]) as f:b=f['amstreea'].arrays(['run','event','tk_pos','tk_dir'],library='np')
        k,found=match({name:v[pick] for name,v in data.items() if name in ['run','event']},b)
        if not found.all():raise ValueError('Geometry identity missing')
        data['tk_pos'][pick]=b['tk_pos'][k];data['tk_dir'][pick]=b['tk_dir'][k]
        source_notes.append(dict(source=int(source),path=lookup[int(source)],events=len(pick)))
        if (num+1)%40==0:print('GEOMETRY_INPUT',num+1,flush=True)
    finite=np.isfinite(data['tk_pos']).all((1,2))&np.isfinite(data['tk_dir']).all((1,2))
    if not finite.all():raise ValueError('Nonfinite native geometry inputs')
    np.savez_compressed(OUT/'measured.npz',**data)
    values=np.column_stack([data['run'],data['event'],data['R'],data['bar'],data['x']*C,data['tk_pos'].reshape(len(ix),-1),data['tk_dir'].reshape(len(ix),-1)])
    formats=['%d','%d','%.12g']+['%d']*4+['%.12g']*(4+54)
    np.savetxt(OUT/'geometry_input.txt',values,fmt=formats)
    protocol=dict(events=len(ix),parent_rigidity_selected=int(ok.sum()),rigidity_GV=[2.5,6],records=source_notes,
        selection='Existing final detector-selected carbon identities intersect native four-used-layer/NewCalib validity domain; deterministic event subsample',
        split='Original whole-run SHA256 modulo5, zero development evaluation, other folds calibration',
        timing='NewCalib LT; native LT is a separate sensitivity branch; time offsets remain fixed',
        RICH_read=False,MC_read=False,absolute_energy_scale_learned=False,
        material='Existing g4mscan map and betaFitter BuildProfile only; no Hybrid calibration configuration loaded',
        hypothesis='Initial model uses two mean-loss scales: upper/lower TOF packages together, and interior material together; not individual active detector deposits',
        approximation='First-order energy-loss timing kernels will be checked against native cell propagation before any deployment claim',
        measured_sha256=digest(OUT/'measured.npz'),geometry_input_sha256=digest(OUT/'geometry_input.txt'))
    (OUT/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');print('PREPARED',len(ix),flush=True)

def moments():
    a=load(OUT/'measured.npz')
    with uproot.open(OUT/'geometry.root') as f:g=f['geometry'].arrays(library='np')
    if not np.array_equal(keys(a),keys(g)):raise ValueError('Native geometry identity mismatch')
    n=len(a['run']);moment=np.zeros((n,4,2,9));column=np.zeros((n,2,9));x=np.zeros((n,4));valid=g['valid'].astype(bool)
    length_error=[]
    for e in range(n):
        if not valid[e]:continue
        hi=np.asarray(g['high'][e],float);lo=np.asarray(g['low'][e],float);length=np.asarray(g['length'][e],float)
        group=np.isin(g['segment'][e],[1,2]).astype(int);ab=np.asarray(g['abundance'][e],float).reshape(-1,9)*ATOMIC_MASS
        nodes=np.r_[hi,lo[-1]];cum=np.r_[0.,np.cumsum(length)]
        center=np.flatnonzero(np.isclose(nodes,g['center_z'][e],rtol=0,atol=1e-5))
        if len(center)!=1:raise ValueError('Reference plane not unique')
        xx=(cum-cum[center[0]])/C
        for l in range(4):
            where=np.flatnonzero(np.isclose(nodes,g['tof_z'][e,l],rtol=0,atol=1e-5))
            if len(where)!=1:raise ValueError('TOF plane not unique')
            x[e,l]=xx[where[0]]
        length_error.append(float(np.max(np.abs(np.diff(x[e])-np.diff(a['x'][e])))))
        middle=.5*(xx[:-1]+xx[1:]);h=np.maximum(x[e,:,None]-middle,0)-np.maximum(-middle,0)-x[e,:,None]*(middle<=0)
        for k in range(2):
            m=group==k;moment[e,:,k]=h[:,m].dot(ab[m]);column[e,k]=ab[m].sum(0)
    if max(length_error)>2e-6:raise ValueError('Native path normalization failed')
    np.savez_compressed(OUT/'moments.npz',run=a['run'],event=a['event'],valid=valid,x=x,moment=moment,column=column)
    result=dict(events=n,valid=int(valid.sum()),max_interlayer_path_over_c_error_ns=max(length_error),
        x_reference='Native L5 plane z=1.698 cm',moment_unit='g cm^-2 ns',column_unit='g cm^-2',
        source='betaFitter HybridEnergyLossCalibration::BuildProfile with native field propagation, voxel boundaries and bar planes',
        group_names=['TOF_packages','inner_material'],geometry_sha256=digest(OUT/'geometry.root'),moments_sha256=digest(OUT/'moments.npz'))
    (OUT/'geometry_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','moments']);p.add_argument('--limit',type=int,default=12000);args=p.parse_args()
    if args.stage=='prepare':prepare(args.limit)
    else:moments()
