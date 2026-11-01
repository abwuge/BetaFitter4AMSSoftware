#!/usr/bin/env python3
"""Repeat the predefined measured-only model on the actual NewCalib LT branch."""
import argparse
import json
import numpy as np
import research_flight_timing_redundancy as r
ROOT=r.ROOT
OLD=r.OUT
OUT=ROOT/'results/flight_calibration_research/newcalib_timing_redundancy'


def prepare():
    OUT.mkdir(exist_ok=True)
    with np.load(str(OLD/'measured.npz'),allow_pickle=False) as f:data={k:f[k] for k in f.files}
    with np.load(str(ROOT/'results/flight_calibration_research/timing_lineage/measured.npz'),allow_pickle=False) as f:
        runs=f['run'];events=f['event'];times=f['tof_tl_NewCalib'][:,0,0,0,0,0,:].astype(float)
    lookup={k:i for i,k in enumerate(zip(runs,events))}
    if len(lookup)!=len(runs):raise ValueError('Duplicate new-calibration identity')
    ix=np.array([lookup[k] for k in zip(data['run'],data['event'])],int)
    if not np.array_equal(runs[ix],data['run']) or not np.array_equal(events[ix],data['event']):raise ValueError('Time branch mismatch')
    new=times[ix];valid=np.isfinite(new).all(axis=1)&(np.abs(new)<1000).all(axis=1)&(new!=0).all(axis=1)
    data['t']=new;data=r.subset(data,valid)
    np.savez_compressed(str(OUT/'measured.npz'),**data)
    r.write(OUT/'protocol.json',dict(baseline='tof_tl_NewCalib[:,0,0,0,0,0,:], LT method0; method1 is identical on audited records',
        baseline_validation='Finite, nonzero, absolute time<1000ns only; this does not establish calibrated quality or absence of pathologies. No new inverse-speed or timing-residual cut.',
        parent_protocol_sha256=r.digest(OLD/'protocol.json'),parent_cache_sha256=r.digest(OLD/'measured.npz'),
        new_time_cache_sha256=r.digest(ROOT/'results/flight_calibration_research/timing_lineage/measured.npz'),
        original_events=len(valid),kept_events=int(valid.sum()),missing_or_numerically_invalid=int((~valid).sum()),
        response_and_correction='Refit independently using the same run cohorts and model family, not copied from native-time coefficients',
        RICH_read=False,MC_coefficients_used=False,comparison='Raw-time redundancy beyond NewCalib LT; no beta resolution claimed'))
    print('NEWCALIB_PREPARED',int(valid.sum()),'events',int((~valid).sum()),'missing/numerically invalid',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','fit','evaluate','bootstrap','all']);p.add_argument('--replicates',type=int,default=32);a=p.parse_args()
    if a.stage in ['prepare','all']:prepare()
    r.OUT=OUT
    if a.stage in ['fit','all']:r.freeze()
    if a.stage in ['evaluate','all']:r.evaluate()
    if a.stage in ['bootstrap','all']:r.bootstrap(a.replicates)
