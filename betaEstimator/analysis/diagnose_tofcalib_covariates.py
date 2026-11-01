#!/usr/bin/env python3
"""Read measured covariates for attribution; never refit the beta model."""
import argparse
import json
import multiprocessing as mp
from pathlib import Path
import numpy as np
import flight_beta_wide as w
from prepare_tofcalib_validation import OUT

FIELDS=['run','event','source','tof_tl','tof_leng','tof_pos','tof_edep','tk_ql',
        'tk_rigidityk','tof_barid','tk_erigidity1','tk_chis1','rich_pb','rich_ebeta1',
        'rich_npe','rich_BetaConsistency','rich_PMTChargeConsistency','rich_pmt',
        'rich_used','rich_pos','rich_q']
LABELS=([f'tof_x{i}' for i in range(4)]+[f'tof_y{i}' for i in range(4)]+
        ['slope_x','slope_y']+[f'path_centered{i}' for i in range(4)]+
        [f'time_residual_over_span{i}' for i in range(4)]+
        [f'log_tof_edep{i}' for i in range(4)]+[f'log_tracker_charge{i}' for i in range(7)]+['log_span'])


def part(source):
    dest=OUT/'parts'/('%05d'%source)
    if not json.loads((dest/'supplement.json').read_text())['selected']:return None
    with np.load(str(dest/'joined.npz'),allow_pickle=False) as f:a={k:f[k] for k in FIELDS}
    w.set_edges(w.R_EDGES);d=w.rt.center_features(a)
    chi=a['tk_chis1'][:,2,2,1,:]
    return dict(run=a['run'],event=a['event'],source=a['source'],continuous=d['continuous'].astype(np.float32),
        bars=a['tof_barid'].astype(np.int16),relR=np.abs(a['tk_erigidity1'][:,2,1]*d['R']),
        chi_x=chi[:,0],chi_y=chi[:,1],rich_pb=a['rich_pb'],rich_error=a['rich_ebeta1'],
        rich_used=a['rich_used'],rich_pmt=a['rich_pmt'],rich_npe_expected=a['rich_npe'][:,2],
        rich_npe_ratio=a['rich_npe'][:,0]/a['rich_npe'][:,2],rich_charge=a['rich_q'][:,1],
        rich_beta_consistency=a['rich_BetaConsistency'],rich_charge_consistency=a['rich_PMTChargeConsistency'],
        rich_x=a['rich_pos'][:,0],rich_y=a['rich_pos'][:,1])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=8);args=parser.parse_args()
    w.check_freeze();indices=json.loads((OUT/'selection.json').read_text())['file_indices'];pieces=[]
    with mp.Pool(args.workers) as pool:
        for i,record in enumerate(pool.imap(part,indices),1):
            if record is not None:pieces.append(record)
            if i%200==0 or i==len(indices):print('Covariates',i,'/',len(indices),flush=True)
    a=w.ev.concatenate(pieces)
    with np.load(str(OUT/'predictions.npz'),allow_pickle=False) as f:
        for k in ('run','event','source'):assert np.array_equal(a[k],f[k]),k
    np.savez_compressed(str(OUT/'diagnostic_covariates.npz'),**a)
    w.ev.dump(OUT/'diagnostic_covariates.json',dict(complete=True,n=len(a['run']),continuous_labels=LABELS,
        predictions_sha256=w.ev.digest(OUT/'predictions.npz'),covariates_sha256=w.ev.digest(OUT/'diagnostic_covariates.npz'),
        script_sha256=w.ev.digest(Path(__file__)),beta_model_refitted=False))


if __name__=='__main__':main()
