#!/usr/bin/env python3
"""Check a real NEWTOF event enriched with newly extracted measurements."""
import json
from pathlib import Path

import numpy as np
import uproot

import flight_beta_wide as w

DEST = w.rt.ROOT / 'results/tofcalib_run_probe'
SOURCE = Path('/lustre02/user/yachen/public/hxwu/newtofcalib_dt_fullsize_mitdst/1305853512_0/1305853512_0.root')
EXTRA = ['tk_rigidityk', 'tk_erigidity1', 'tof_bpatt', 'rich_itrtrack', 'rich_ebeta1']


def ids(data):
    return np.rec.fromarrays([data['run'].astype(np.int64), data['event'].astype(np.int64)],
                            names='run,event')


def main():
    w.check_freeze()
    names = list(dict.fromkeys(w.ev.MEASURED + w.ev.SELECTION + ['rich_q']))
    with uproot.open(str(DEST / 'carbon_event.root')) as f:
        native = f['amstreea'].arrays(names, library='np')
    if len(native['run']) != 1 or int(native['run'][0]) != 1305853512 or int(native['event'][0]) != 14410:
        raise RuntimeError('Target event was not produced')
    with uproot.open(str(SOURCE)) as f:
        original = f['amstreea'].arrays([k for k in names if k not in EXTRA], library='np')
    common, original_rows, extra_rows = np.intersect1d(ids(original), ids(native), return_indices=True)
    if len(common) != len(native['run']):
        raise RuntimeError('Incomplete event join')
    joined = w.ev.sub(original, original_rows)
    for key in EXTRA:
        joined[key] = native[key][extra_rows]
    # Persist only the genuinely missing measurements in a sidecar. Do not
    # replace the NEWTOF times/geometry with the probe producer's output.
    supplement = {k: native[k][extra_rows] for k in ['run', 'event'] + EXTRA}
    with uproot.recreate(str(DEST / 'missing_branches.root')) as f:
        f['beta_extra'] = supplement
    with uproot.open(str(DEST / 'missing_branches.root')) as f:
        replay = f['beta_extra'].arrays(library='np')
    for key in supplement:
        if not np.array_equal(replay[key], supplement[key]):
            raise RuntimeError('Supplement serialization changed a measurement')
    w.set_edges(w.R_EDGES)
    selected = w.selection(joined)
    ring = w.rt.rich_quality(joined) & (joined['rich_q'][:, 1] > 4) & (joined['rich_q'][:, 1] < 8)
    if not selected.all() or not ring.all():
        raise RuntimeError('Enriched NEWTOF event failed current selection')
    data = w.rt.center_features(joined)
    corrected = w.predict(data)
    if not np.isfinite(corrected).all():
        raise RuntimeError('Frozen model could not predict enriched NEWTOF')
    geometry = w.rt.MaterialMap().build_columns(joined['tk_pos'], joined['tk_dir'])
    mass = 11.17793
    expected_u = np.sqrt(1 + (mass / (6 * data['R'])) ** 2)
    physical = w.rt.propagate(expected_u, data['R'], geometry, mass_override=mass)
    delay = 1 / physical['beta_planes'][:, -1] - expected_u
    if not physical['valid'].all():
        raise RuntimeError('Material propagation failed')
    old, new = 1 / (data['u'] + delay), 1 / (corrected + delay)
    same_times = bool(np.array_equal(joined['tof_tl'], original['tof_tl'][original_rows]))
    if not same_times:
        raise RuntimeError('NEWTOF time was overwritten')
    record = dict(source=str(SOURCE), target_run=1305853512, target_event=14410,
                  missing_branch_keys=EXTRA, join_complete=True, selected=True, rich_reference=True,
                  original_newtof_time_preserved=same_times, time_branch='tof_tl',
                  native_time_equal=bool(np.array_equal(joined['tof_tl'], native['tof_tl'][extra_rows])),
                  corrected_u=corrected, measured_rigidity_GV=data['R'],
                  baseline_beta_rich=old, corrected_beta_rich=new, measured_beta_rich=joined['rich_beta'][:, 0],
                  propagation_valid=True, training_or_recalibration_performed=False,
                  performance_claim='one-event compatibility probe; no statistical performance inference',
                  model_freeze_sha256=w.ev.digest(w.DEST / 'freeze.json'))
    w.ev.dump(DEST / 'joined_verification.json', record)
    print('PASS: actual NEWTOF + five measured supplement branches -> selection -> frozen model -> RICH evaluation', flush=True)
    print('run/event',record['target_run'],record['target_event'],'R',record['measured_rigidity_GV'],flush=True)


if __name__ == '__main__':
    main()
