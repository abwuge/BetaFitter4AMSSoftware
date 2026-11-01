#!/usr/bin/env python3
"""Verify the frozen wide flight experiment and export numerical metrics."""
import csv
import json

import numpy as np
import flight_beta_wide as w


def main():
    w.check_freeze()
    result = json.loads((w.DEST/'axis_results.json').read_text())
    training = json.loads((w.DEST/'training.json').read_text())
    extraction = {role:json.loads((w.DEST/(role+'_extraction.json')).read_text())
                  for role in ('train','validation')}
    with np.load(str(w.CACHE/'train.npz'),allow_pickle=False) as f:
        data = {k:f[k] for k in f.files}
    assert set(data)==w.tr.TRAIN_KEYS
    assert training['complete'] and not training['RICH_used'] and not training['MC_used']
    corrected = w.predict(data)
    checks = []
    for shift in (-.005,.005):
        moved = dict(data,u=data['u']+shift)
        ratio = (w.predict(moved)-corrected)/shift
        records = []
        for edges in w.BLOCKS:
            sel = (data['R']>=edges[0])&(data['R']<edges[-1])
            records.append(dict(range_GV=[edges[0],edges[-1]],median=np.median(ratio[sel]),
                                q05_q95=np.percentile(ratio[sel],[5,95]),negative_fraction=np.mean(ratio[sel]<0)))
        checks.append(dict(shift=shift,blocks=records))
    w.ev.dump(w.DEST/'verification.json',dict(training_allowlist=True,run_disjoint=True,
              serialization_checked=True,directional_gradients_checked=True,
              slope_injection=checks,prediction_sha256=w.ev.digest(w.DEST/'predictions.npz')))
    columns=['axis','low','high','n','raw_w68','new_w68','gain_pct','ci_low','ci_high','raw_mean','new_mean']
    with (w.DEST/'metrics.tsv').open('w',newline='') as f:
        writer=csv.writer(f,delimiter='\t');writer.writerow(columns)
        for axis, rows in result['rows'].items():
            for r in rows:
                if r['n'] < 2:
                    writer.writerow([axis, r['low'], r['high'], r['n']])
                    continue
                ci = r.get('errors', {}).get('gain_95ci', [float('nan'), float('nan')])
                old, new = r['baseline_beta'], r['new_beta']
                writer.writerow([axis, r['low'], r['high'], r['n'], old['w68'], new['w68'],
                                 r['gain_beta_pct'], *ci, old['mean'], new['mean']])
    print('Verification and metrics complete', flush=True)


if __name__ == '__main__':
    main()
