#!/usr/bin/env python3
"""Join prospective measured supplements and build the frozen input domain."""
import argparse
import json
from pathlib import Path
import numpy as np
import uproot
from prepare_prospective_flight import ROOT, OUT
from research_flight_conditional import load, write, digest
from research_conditional_mc import stamp
import research_flight_material as material


def plan():
    path = OUT / 'native_domain_protocol.json'
    if path.exists():
        raise FileExistsError('Native measurement-domain protocol already exists')
    write(path, dict(created_utc=stamp(), source_sha256=digest(Path(__file__)),
        candidates_sha256=digest(OUT / 'candidates.npz'), preparation_sha256=digest(OUT / 'preparation.json'),
        source_domain='Previously fixed carbon, central Tracker trajectory, positive ADC and original four-used finite native/NewCalib layer-time domain',
        native_domain='Native Track/BetaH association matches; GetSumHit=GetUseHit=4, all four used and good times, BetaPattern=4444; 2.5<=R<6 GV, 0<abs(GetErrRinv*R)<0.2; both native normalized Tracker chi-square X/Y in (0,10]',
        R='GetRigidity(fit,1), the same z=0 Kalman convention used by the frozen models; the separately interpolated L5 value is audit-only',
        timing='Stored native/NewCalib times remain authoritative. Native time, associated end ADC and bar must match on every native selected event before applying R or chi-square cuts.',
        paths='Use stored layer path; require native interlayer paths within 1e-5 cm on every native selected event',
        unmatched='Report every native rejection; do not alter association to improve beta agreement',
        geometry='Existing native BuildProfile and unchanged material moment exporter. Do not load old fitted Hybrid coefficients.',
        geometry_source_sha256=digest(ROOT / 'analysis/flight_material_geometry.C'),
        geometry_binary_sha256=digest(ROOT / 'results/flight_calibration_research/mean_loss_likelihood/supplement'),
        moment_source_sha256=digest(ROOT / 'analysis/research_flight_material.py'),
        no_output_beta_cut=True, new_RICH_values_read=False))
    print('PROSPECTIVE_NATIVE_DOMAIN_PROTOCOL', flush=True)


def check():
    p = json.loads((OUT / 'native_domain_protocol.json').read_text())
    for path, key in [(Path(__file__), 'source_sha256'), (OUT / 'candidates.npz', 'candidates_sha256'),
                      (OUT / 'preparation.json', 'preparation_sha256'), (ROOT / 'analysis/research_flight_material.py', 'moment_source_sha256')]:
        if digest(path) != p[key]:
            raise ValueError('Prospective measured pipeline changed: ' + key)
    return p


def join():
    check()
    if (OUT / 'measured.npz').exists():
        raise FileExistsError('Prospective measured cache already exists')
    a = load(OUT / 'candidates.npz'); prep = json.loads((OUT / 'preparation.json').read_text())
    batch = json.loads((OUT / 'batch/protocol.json').read_text()); n = len(a['run'])
    seen = np.zeros(n, bool); selected = np.zeros(n, bool); outputs = {}; records = []
    for row in prep['jobs']:
        source = row['source']; ix = np.flatnonzero(a['source'] == source); folder = OUT / 'native' / ('%05d' % source)
        done = json.loads((folder / 'complete.json').read_text())
        if not done['complete'] or done['binary_sha256'] != batch['binary_sha256'] or digest(folder / 'measurement.root') != done['output_sha256']:
            raise ValueError('Incomplete or changed native measurement output')
        with uproot.open(folder / 'measurement.root') as f:
            native = f['measured'].arrays(library='np')
        if len(ix) != len(native['run']) or not native['found'].all():
            raise ValueError('Missing native candidate identity')
        for key in ['run', 'event', 'itrtrack', 'requested_beta']:
            original = a['ibetah'][ix] if key == 'requested_beta' else a[key][ix]
            if not np.array_equal(native[key], original):
                raise ValueError('Native candidate identity or requested association differs')
        selected[ix] = native['selected'].astype(bool)
        for key, value in native.items():
            if key not in outputs:
                outputs[key] = np.empty((n,) + value.shape[1:], value.dtype)
            outputs[key][ix] = value
        seen[ix] = True
        records.append(dict(source=source, events=len(ix), native_selected=int(native['selected'].sum()),
                            completion_sha256=digest(folder / 'complete.json'), output_sha256=done['output_sha256']))
    if not seen.all():
        raise ValueError('Prospective candidate join incomplete')
    b = outputs; comparisons = {}
    if not np.array_equal(a['ibetah'][selected], b['ibetah'][selected]):
        raise ValueError('Selected native BetaH association differs')
    for name, old, new in [('time', a['t_native'], b['time']), ('ADC', a['adc'], b['adc']), ('bar', a['bar'], b['bar']),
                           ('path_cm', a['x'] * 29.9792458, b['length']), ('chi', a['tk_chis1'][:, 2, 2, 1, :], b['chi'])]:
        diff = np.asarray(new[selected], float) - np.asarray(old[selected], float)
        comparisons[name] = dict(max_abs=float(np.max(abs(diff))), exact=bool(np.array_equal(old[selected], new[selected])),
                                 rms=float(np.sqrt(np.mean(diff * diff))))
    if not all(comparisons[k]['exact'] for k in ['time', 'ADC', 'bar']) or comparisons['path_cm']['max_abs'] > 1e-5:
        raise ValueError('Native readout equivalence failed: ' + str(comparisons))
    R = b['R'][:, 1].astype(float); sr = np.abs(b['sigma_r'].astype(float)); relative = np.abs(sr * R)
    keep = selected.copy(); counts = [dict(stage='native_association_and_fit', events=int(keep.sum()))]
    conditions = [('four_used_good_pattern', (b['sumhit'] == 4) & (b['usehit'] == 4) & b['used'].all(1) & b['good'].all(1) & (b['bpatt'] == 4444)),
                  ('rigidity_2p5_6GV', np.isfinite(R) & (R >= 2.5) & (R < 6)),
                  ('relative_R_error', np.isfinite(relative) & (relative > 0) & (relative < .2)),
                  ('tracker_chi', np.isfinite(b['chi'][:, :2]).all(1) & ((b['chi'][:, :2] > 0) & (b['chi'][:, :2] <= 10)).all(1))]
    for name, mask in conditions:
        keep &= mask; counts.append(dict(stage=name, events=int(keep.sum())))
    if not keep.any():
        raise ValueError('Empty final prospective native input domain')
    with np.load(OUT / 'excluded_runs.npz') as old:
        overlap = np.intersect1d(a['run'][keep], old['run'])
    if len(overlap):
        raise ValueError('Historical run entered prospective detector domain')
    native = {k: v[keep] for k, v in a.items()}
    native.update(R=R[keep], relR=relative[keep], sigma_r=sr[keep], R_L5=b['R_L5'][keep],
                  tof_betah=b['beta'][keep], native_chi=b['chi'][keep], native_time_error=b['etime'][keep],
                  good=b['good'][keep].astype(bool))
    np.savez_compressed(OUT / 'measured.npz', **native)
    values = np.column_stack([native['run'], native['event'], native['R'], native['bar'], native['x'] * 29.9792458,
                              native['tk_pos'].reshape(keep.sum(), -1), native['tk_dir'].reshape(keep.sum(), -1)])
    np.savetxt(OUT / 'geometry_input.txt', values, fmt=['%d', '%d', '%.12g'] + ['%d'] * 4 + ['%.12g'] * 58)
    np.savez_compressed(OUT / 'native_selection.npz', run=a['run'], event=a['event'], selected=selected, final_domain=keep)
    plane = (native['R_L5'] - native['R']) / native['R']
    write(OUT / 'native_verification.json', dict(complete=True, created_utc=stamp(), candidates=n, events=int(keep.sum()),
        runs=len(np.unique(native['run'])), historical_run_overlap=len(overlap), stage_counts=counts, comparisons=comparisons, records=records,
        measured_sha256=digest(OUT / 'measured.npz'), geometry_input_sha256=digest(OUT / 'geometry_input.txt'),
        selection_sha256=digest(OUT / 'native_selection.npz'), native_protocol_sha256=digest(OUT / 'native_domain_protocol.json'),
        relative_L5_minus_z0_R_percentiles=np.percentile(plane, [0, 5, 50, 95, 100]),
        original_measurements_preserved=True, new_RICH_values_read=False))
    print('PROSPECTIVE_NATIVE_JOINED', int(keep.sum()), len(np.unique(native['run'])), counts, comparisons, flush=True)


def moments():
    check()
    if (OUT / 'moments.npz').exists():
        raise FileExistsError('Prospective material moments already exist')
    material.OUT = OUT
    material.moments()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('stage', choices=['plan', 'join', 'moments']); args = parser.parse_args()
    globals()[args.stage]()
