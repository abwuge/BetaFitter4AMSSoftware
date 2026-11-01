#!/usr/bin/env python3
"""Open the new RICH reference only after the prospective predictions are frozen."""
import argparse
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import uproot
from prepare_prospective_flight import ROOT, OUT
from predict_prospective_flight import check as prediction_check
from research_flight_conditional import load, write, digest
from research_conditional_mc import stamp
import rich_transport as material

REF = OUT / 'reference'
BINARY = OUT / 'bin/rich_reference'
FIELDS = ['run', 'event', 'itrtrack', 'irich', 'rich_itrtrack', 'rich_beta',
          'rich_good', 'rich_clean', 'rich_pb', 'rich_pmt', 'rich_used',
          'rich_BetaConsistency', 'rich_PMTChargeConsistency', 'rich_npe',
          'rich_NaF', 'rich_pos', 'rich_tile']


def gate():
    prediction_check()
    frozen = json.loads((OUT / 'output_freeze.json').read_text())
    if not frozen['complete'] or frozen['new_RICH_values_read']:
        raise ValueError('A complete reference-blind prediction snapshot is required')
    for name, key in [('predictions_before_reference.npz', 'prediction_sha256'),
                      ('measured.npz', 'measured_sha256'), ('moments.npz', 'moments_sha256'),
                      ('rich_geometry.npz', 'rich_geometry_sha256'), ('prediction_protocol.json', 'protocol_sha256')]:
        if digest(OUT / name) != frozen[key]:
            raise ValueError('Prospective frozen input changed: ' + name)
    return frozen


def plan():
    frozen = gate()
    REF.mkdir(exist_ok=True)
    if (REF / 'protocol.json').exists():
        raise FileExistsError('Prospective reference protocol already exists')
    write(REF / 'protocol.json', dict(created_utc=stamp(), events=frozen['events'],
        output_freeze_sha256=digest(OUT / 'output_freeze.json'), source_sha256=digest(Path(__file__)),
        quality_source_sha256=digest(Path(material.__file__)), native_binary_sha256=digest(BINARY),
        native_source_sha256=digest(ROOT / 'analysis/prospective_rich_reference.C'),
        ntuple_producer_sha256=digest(ROOT.parent / 'amsd69n/analysis_amsd69n.C'),
        npe_columns=['ring collected photoelectrons', 'ring expected photoelectrons', 'all-event collected photoelectrons'],
        reference_fields=FIELDS, optional_field='rich_ebeta1 is recorded when present; the independently reread native getBetaError is authoritative',
        reference='Native existing associated RichRing getBeta(), checked against stored rich_beta[0] within 1e-7 before any residual is evaluated',
        association='Exact run/event and requested Tracker/ring index; ring iTrTrack must match the fixed measured Tracker',
        quality='Unchanged rich_transport.rich_quality(strict=True), NaF, plus finite native total Tracker chi-square in (0,10]',
        cuts=dict(beta=[.78, 1.02], beta_error=[0., .003], probability_min=.05, PMTs_min=4, used_hits_min=6,
                  beta_consistency_max_NaF=.01, PMT_charge_consistency_max=15., ring_to_event_collected_npe_min=.4,
                  expected_npe_min_NaF=1., NaF_max_abs_xy_cm=17.),
        native_preselection='Apply all the fixed parent-reference quality cuts except beta error, and fixed measured total Tracker chi-square; native reread only this prespecified parent-qualified cohort',
        missing_reference='Keep all measured identities and reference-valid flags; ring absence or bad RICH quality is a reference rejection, not a reconstructed-output rejection',
        exposure='Conservatively record every source file whose RICH branches are opened; independence is relative to the persisted histories enumerated in exposure.json',
        no_TOF_RICH_agreement_cut=True, no_reconstructed_beta_cut=True, no_mass_cut=True,
        no_model_or_threshold_refit=True, new_RICH_values_read=False))
    print('PROSPECTIVE_REFERENCE_PROTOCOL', frozen['events'], flush=True)


def check():
    gate()
    p = json.loads((REF / 'protocol.json').read_text())
    for path, key in [(Path(__file__), 'source_sha256'), (Path(material.__file__), 'quality_source_sha256'),
                      (BINARY, 'native_binary_sha256'), (OUT / 'output_freeze.json', 'output_freeze_sha256')]:
        if digest(path) != p[key]:
            raise ValueError('Prospective reference implementation changed: ' + key)
    return p


def prepare():
    p = check()
    if (REF / 'parent.npz').exists() or (REF / 'preparation.json').exists():
        raise FileExistsError('New RICH parent values were already prepared')
    if (REF / 'read_started.json').exists():
        raise FileExistsError('Inspect the previously started reference read before retrying')
    a = load(OUT / 'measured.npz')
    original = json.loads((OUT / 'preparation.json').read_text())
    source_map = {r['source']: r for r in original['records']}
    n = len(a['run'])
    sources = sorted(int(s) for s in np.unique(a['source']))
    write(REF / 'read_started.json', dict(created_utc=stamp(), protocol_sha256=digest(REF / 'protocol.json'),
        output_freeze_sha256=p['output_freeze_sha256'], files=[source_map[s]['path'] for s in sources],
        stage='First opening of prospective RICH branch values; model predictions were already frozen'))

    def one(source):
        ix = np.flatnonzero(a['source'] == source)
        entries = a['entry'][ix]
        start, stop = int(entries.min()), int(entries.max()) + 1
        with uproot.open(source_map[source]['path']) as f:
            tree = f['amstreea']
            optional = ['rich_ebeta1'] if 'rich_ebeta1' in tree else []
            raw = tree.arrays(FIELDS + optional, entry_start=start, entry_stop=stop, library='np')
        out = {k: v[entries - start] for k, v in raw.items()}
        for key in ['run', 'event', 'itrtrack']:
            if not np.array_equal(out[key], a[key][ix]):
                raise ValueError('New RICH parent identity differs: ' + key)
        if 'rich_ebeta1' not in out:
            out['rich_ebeta1'] = np.full(len(ix), np.nan)
        return ix, out, dict(source=source, path=source_map[source]['path'], events=len(ix),
                             opened_entry_range=[start, stop], stored_error_present=bool(optional))

    parent = {}
    seen = np.zeros(n, bool)
    evidence = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        for ix, row, note in pool.map(one, sources):
            for key, value in row.items():
                if key not in parent:
                    parent[key] = np.empty((n,) + value.shape[1:], value.dtype)
                parent[key][ix] = value
            seen[ix] = True
            evidence.append(note)
    if not seen.all():
        raise ValueError('Incomplete parent-reference join')
    trial = dict(parent)
    trial['rich_ebeta1'] = np.full(n, .001)
    chi = a['native_chi'][:, 2]
    measured_quality = np.isfinite(chi) & (chi > 0) & (chi <= 10)
    pre = material.rich_quality(trial, strict=True) & parent['rich_NaF'].astype(bool) & measured_quality
    parent.update(source=a['source'], entry=a['entry'], pre_native=pre, measured_total_chi_quality=measured_quality)
    np.savez_compressed(REF / 'parent.npz', **parent)
    (REF / 'targets').mkdir(exist_ok=True)
    (REF / 'parts').mkdir(exist_ok=True)
    (REF / 'claims').mkdir(exist_ok=True)
    jobs = []
    for source in sources:
        ix = np.flatnonzero((a['source'] == source) & pre)
        if not len(ix):
            continue
        target = REF / 'targets' / ('%05d.txt' % source)
        np.savetxt(target, np.column_stack([a['run'][ix], a['event'][ix], a['itrtrack'][ix], parent['irich'][ix]]), fmt='%d')
        raw = OUT / 'native' / ('%05d' % source) / 'raw_files.txt'
        jobs.append(dict(source=source, events=len(ix), targets=str(target), targets_sha256=digest(target),
                         raw_files=str(raw), raw_files_sha256=digest(raw)))
    write(REF / 'preparation.json', dict(complete=True, created_utc=stamp(), measured_events=n,
        parent_qualified_events=int(pre.sum()), jobs=jobs, opened_sources=evidence,
        parent_sha256=digest(REF / 'parent.npz'), protocol_sha256=digest(REF / 'protocol.json'),
        new_RICH_values_read=True, predictions_changed=False))
    print('PROSPECTIVE_RICH_PARENT', n, 'native_targets', int(pre.sum()), 'sources', len(jobs), flush=True)


def collect():
    p = check()
    preparation = json.loads((REF / 'preparation.json').read_text())
    if (REF / 'reference.npz').exists():
        raise FileExistsError('Prospective RICH reference already collected')
    if digest(REF / 'parent.npz') != preparation['parent_sha256']:
        raise ValueError('Parent-reference cache changed')
    parent = load(REF / 'parent.npz')
    n = len(parent['run'])
    seen = np.zeros(n, bool)
    found = np.zeros(n, bool)
    associated = np.zeros(n, bool)
    beta = np.full((n, 3), np.nan)
    error = np.full(n, np.nan)
    parts = {}
    differences = []
    for row in preparation['jobs']:
        path = REF / 'parts' / ('%05d.root' % row['source'])
        note = json.loads(path.with_suffix('.json').read_text())
        if not note['complete'] or digest(path) != note['output_sha256'] or note['binary_sha256'] != p['native_binary_sha256']:
            raise ValueError('Incomplete or changed native RICH part')
        with uproot.open(path) as f:
            native = f['reference'].arrays(library='np')
        ix = np.flatnonzero((parent['source'] == row['source']) & parent['pre_native'])
        for key in ['run', 'event', 'itrtrack', 'irich']:
            if not np.array_equal(native[key], parent[key][ix]):
                raise ValueError('Native RICH target identity mismatch: ' + key)
        if not native['found'].all():
            raise ValueError('Missing raw event in new RICH native reread')
        seen[ix] = True
        found[ix] = native['ring_found'].astype(bool)
        associated[ix] = native['rich_itrtrack'] == parent['itrtrack'][ix]
        beta[ix] = native['beta']
        error[ix] = native['error']
        same = found[ix] & associated[ix]
        diff = native['beta'][same, 0].astype(float) - parent['rich_beta'][ix[same], 0].astype(float)
        if not np.isfinite(diff).all() or np.any(abs(diff) > 1e-7):
            raise ValueError('Stored/native RICH beta mismatch; audit it before scoring')
        differences.extend(diff.tolist())
        parts[str(row['source'])] = dict(output=note['output_sha256'], completion=digest(path.with_suffix('.json')))
    if not np.array_equal(seen, parent['pre_native']):
        raise ValueError('Native RICH reread does not cover the fixed parent-qualified cohort')
    checked = dict(parent)
    checked['rich_beta'] = beta
    checked['rich_ebeta1'] = error
    valid = seen & found & associated & parent['measured_total_chi_quality'] & parent['rich_NaF'].astype(bool)
    valid &= material.rich_quality(checked, strict=True)
    reason = np.full(n, 'parent_quality_or_measured_chi', dtype='U40')
    reason[seen] = 'native_RICH_quality'
    reason[seen & ~found] = 'native_ring_missing'
    reason[seen & found & ~associated] = 'native_ring_track_mismatch'
    reason[valid] = 'accepted'
    np.savez_compressed(REF / 'reference.npz', run=parent['run'], event=parent['event'], source=parent['source'],
        entry=parent['entry'], valid=valid, beta=beta[:, 0], beta_variants=beta, error=error, reason=reason)
    delta = np.asarray(differences)
    write(REF / 'audit.json', dict(complete=True, created_utc=stamp(), measured_events=n,
        native_requested=int(seen.sum()), accepted=int(valid.sum()), runs=len(np.unique(parent['run'][valid])),
        rejection_counts={str(k): int(v) for k, v in zip(*np.unique(reason, return_counts=True))},
        native_stored_beta_max_abs=float(np.max(abs(delta))) if len(delta) else None,
        native_stored_beta_exact=bool(np.all(delta == 0)),
        accepted_error_percentiles=np.percentile(error[valid], [0, 5, 50, 95, 100]) if valid.any() else [],
        reference_sha256=digest(REF / 'reference.npz'), protocol_sha256=digest(REF / 'protocol.json'),
        preparation_sha256=digest(REF / 'preparation.json'), output_freeze_sha256=p['output_freeze_sha256'],
        parts=parts, new_RICH_values_read=True, no_output_based_selection=True))
    print('PROSPECTIVE_RICH_REFERENCE', int(valid.sum()), 'of', n, 'runs', len(np.unique(parent['run'][valid])), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['plan', 'prepare', 'collect'])
    globals()[parser.parse_args().stage]()
