#!/usr/bin/env python3
"""Native association completion for NEWTOF trees lacking rich_itrtrack.

The original prediction, score and reference protocols remain immutable. This
reader records a schema amendment and obtains the missing association from the
already frozen native ring reader, rather than pretending an ntuple value exists.
"""
import argparse
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import uproot
import reference_prospective_flight as original
from reference_prospective_flight import ROOT, OUT, REF, BINARY
from research_flight_conditional import load, write, digest
from research_conditional_mc import stamp
import rich_transport as material

FIELDS = [key for key in original.FIELDS if key != 'rich_itrtrack']


def amend():
    p = original.check()
    if (REF / 'schema_amendment.json').exists():
        raise FileExistsError('Reference schema amendment already recorded')
    if (REF / 'parent.npz').exists() or (REF / 'preparation.json').exists():
        raise ValueError('Do not amend an already prepared reference population')
    if not (REF / 'read_started.json').exists():
        raise ValueError('This amendment documents the actual original schema failure')
    write(REF / 'schema_amendment.json', dict(created_utc=stamp(), complete=True,
        reason='Original uproot read failed because the NEWTOF ntuple has no rich_itrtrack branch; its absence was verified from branch metadata, not a beta residual.',
        original_protocol_sha256=digest(REF / 'protocol.json'), original_reader_sha256=digest(Path(original.__file__)),
        source_sha256=digest(Path(__file__)), original_attempt_sha256=digest(REF / 'read_started.json'),
        native_binary_sha256=p['native_binary_sha256'], prediction_freeze_sha256=digest(OUT / 'output_freeze.json'),
        correction='Defer an absent stored ring-to-Tracker association to the native RichRing iTrTrack reread. Unknown stored values are explicitly marked unavailable, never treated as measured matches.',
        invariant='Identical requested ring and Tracker, identical numerical RICH quality cuts, identical native getBetaError requirement, identical model predictions, common-event definition and comparison endpoints.',
        parent_preselection='All existing parent quality predicates except genuinely absent association/error fields; the final quality decision uses observed native association and error.',
        exposure='Retain the original full-source exposure ledger even though the first read failed during expression resolution. No claim that the attempted sources are untouched.',
        beta_residuals_evaluated_before_amendment=False, model_or_threshold_change=False))
    print('PROSPECTIVE_REFERENCE_SCHEMA_AMENDMENT', flush=True)


def check():
    p = original.check()
    amendment = json.loads((REF / 'schema_amendment.json').read_text())
    for path, key in [(Path(__file__), 'source_sha256'), (Path(original.__file__), 'original_reader_sha256'),
                      (REF / 'protocol.json', 'original_protocol_sha256'), (REF / 'read_started.json', 'original_attempt_sha256'),
                      (BINARY, 'native_binary_sha256'), (OUT / 'output_freeze.json', 'prediction_freeze_sha256')]:
        if digest(path) != amendment[key]:
            raise ValueError('Reference schema amendment input changed: ' + key)
    return p


def prepare():
    p = check()
    if (REF / 'parent.npz').exists() or (REF / 'preparation.json').exists() or (REF / 'schema_read_started.json').exists():
        raise FileExistsError('Inspect the previously started amended reference read before retrying')
    a = load(OUT / 'measured.npz')
    preparation = json.loads((OUT / 'preparation.json').read_text())
    source_map = {r['source']: r for r in preparation['records']}
    n = len(a['run'])
    sources = sorted(int(s) for s in np.unique(a['source']))
    write(REF / 'schema_read_started.json', dict(created_utc=stamp(), schema_amendment_sha256=digest(REF / 'schema_amendment.json'),
        output_freeze_sha256=p['output_freeze_sha256'], sources=[source_map[s]['path'] for s in sources]))

    def one(source):
        ix = np.flatnonzero(a['source'] == source)
        entries = a['entry'][ix]
        start, stop = int(entries.min()), int(entries.max()) + 1
        with uproot.open(source_map[source]['path']) as f:
            tree = f['amstreea']
            optional = [key for key in ['rich_ebeta1', 'rich_itrtrack'] if key in tree]
            raw = tree.arrays(FIELDS + optional, entry_start=start, entry_stop=stop, library='np')
        out = {key: value[entries - start] for key, value in raw.items()}
        for key in ['run', 'event', 'itrtrack']:
            if not np.array_equal(out[key], a[key][ix]):
                raise ValueError('New RICH parent identity differs: ' + key)
        out['stored_association_available'] = np.full(len(ix), 'rich_itrtrack' in optional, bool)
        if 'rich_itrtrack' not in out:
            out['rich_itrtrack'] = np.full(len(ix), -1, np.int32)
        if 'rich_ebeta1' not in out:
            out['rich_ebeta1'] = np.full(len(ix), np.nan)
        return ix, out, dict(source=source, path=source_map[source]['path'], events=len(ix),
            opened_entry_range=[start, stop], stored_error_present='rich_ebeta1' in optional,
            stored_association_present='rich_itrtrack' in optional)

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
        raise ValueError('Incomplete amended parent-reference join')
    # These temporary replacements defer unavailable predicates only. They are
    # not saved as measured ring associations or measured error estimates.
    trial = dict(parent)
    trial['rich_ebeta1'] = np.full(n, .001)
    trial['rich_itrtrack'] = np.where(parent['stored_association_available'], parent['rich_itrtrack'], parent['itrtrack'])
    chi = a['native_chi'][:, 2]
    measured_quality = np.isfinite(chi) & (chi > 0) & (chi <= 10)
    pre = material.rich_quality(trial, strict=True) & parent['rich_NaF'].astype(bool) & measured_quality
    parent.update(source=a['source'], entry=a['entry'], pre_native=pre, measured_total_chi_quality=measured_quality)
    np.savez_compressed(REF / 'parent.npz', **parent)
    for name in ['targets', 'parts', 'claims']:
        (REF / name).mkdir(exist_ok=True)
    jobs = []
    for source in sources:
        ix = np.flatnonzero((a['source'] == source) & pre)
        if not len(ix):
            continue
        path = REF / 'targets' / ('%05d.txt' % source)
        np.savetxt(path, np.column_stack([a['run'][ix], a['event'][ix], a['itrtrack'][ix], parent['irich'][ix]]), fmt='%d')
        raw = OUT / 'native' / ('%05d' % source) / 'raw_files.txt'
        jobs.append(dict(source=source, events=len(ix), targets=str(path), targets_sha256=digest(path),
            raw_files=str(raw), raw_files_sha256=digest(raw)))
    write(REF / 'preparation.json', dict(complete=True, created_utc=stamp(), measured_events=n,
        parent_qualified_events=int(pre.sum()), jobs=jobs, opened_sources=evidence,
        missing_stored_associations=int((~parent['stored_association_available']).sum()),
        parent_sha256=digest(REF / 'parent.npz'), protocol_sha256=digest(REF / 'protocol.json'),
        schema_amendment_sha256=digest(REF / 'schema_amendment.json'),
        new_RICH_values_read=True, predictions_changed=False))
    print('PROSPECTIVE_RICH_PARENT', n, 'native_targets', int(pre.sum()), 'sources', len(jobs), flush=True)


def collect():
    p = check()
    preparation = json.loads((REF / 'preparation.json').read_text())
    if (REF / 'reference.npz').exists():
        raise FileExistsError('Prospective reference already collected')
    if digest(REF / 'parent.npz') != preparation['parent_sha256'] or digest(REF / 'schema_amendment.json') != preparation['schema_amendment_sha256']:
        raise ValueError('Amended reference preparation changed')
    parent = load(REF / 'parent.npz')
    n = len(parent['run'])
    seen = np.zeros(n, bool)
    found = np.zeros(n, bool)
    association = np.full(n, -1, np.int32)
    beta = np.full((n, 3), np.nan)
    error = np.full(n, np.nan)
    parts = {}
    differences = []
    for row in preparation['jobs']:
        path = REF / 'parts' / ('%05d.root' % row['source'])
        note = json.loads(path.with_suffix('.json').read_text())
        if not note['complete'] or digest(path) != note['output_sha256'] or note['binary_sha256'] != p['native_binary_sha256']:
            raise ValueError('Incomplete or changed native reference part')
        with uproot.open(path) as f:
            native = f['reference'].arrays(library='np')
        ix = np.flatnonzero((parent['source'] == row['source']) & parent['pre_native'])
        for key in ['run', 'event', 'itrtrack', 'irich']:
            if not np.array_equal(native[key], parent[key][ix]):
                raise ValueError('Native RICH target identity mismatch: ' + key)
        if not native['found'].all():
            raise ValueError('Raw target event missing from native RICH reread')
        seen[ix] = True
        found[ix] = native['ring_found'].astype(bool)
        association[ix] = native['rich_itrtrack']
        beta[ix] = native['beta']
        error[ix] = native['error']
        same = found[ix] & (association[ix] == parent['itrtrack'][ix])
        diff = native['beta'][same, 0].astype(float) - parent['rich_beta'][ix[same], 0].astype(float)
        if not np.isfinite(diff).all() or np.any(abs(diff) > 1e-7):
            raise ValueError('Stored/native RICH beta mismatch; audit it before computing model residuals')
        differences.extend(diff.tolist())
        parts[str(row['source'])] = dict(output=note['output_sha256'], completion=digest(path.with_suffix('.json')))
    if not np.array_equal(seen, parent['pre_native']):
        raise ValueError('Native RICH reread does not cover the fixed parent-qualified cohort')
    associated = association == parent['itrtrack']
    checked = dict(parent)
    checked['rich_beta'] = beta
    checked['rich_ebeta1'] = error
    checked['rich_itrtrack'] = association
    valid = seen & found & associated & parent['measured_total_chi_quality'] & parent['rich_NaF'].astype(bool)
    valid &= material.rich_quality(checked, strict=True)
    reason = np.full(n, 'parent_quality_or_measured_chi', dtype='U40')
    reason[seen] = 'native_RICH_quality'
    reason[seen & ~found] = 'native_ring_missing'
    reason[seen & found & ~associated] = 'native_ring_track_mismatch'
    reason[valid] = 'accepted'
    np.savez_compressed(REF / 'reference.npz', run=parent['run'], event=parent['event'], source=parent['source'], entry=parent['entry'],
        valid=valid, beta=beta[:, 0], beta_variants=beta, error=error, native_rich_itrtrack=association, reason=reason)
    diff = np.asarray(differences)
    write(REF / 'audit.json', dict(complete=True, created_utc=stamp(), measured_events=n,
        native_requested=int(seen.sum()), accepted=int(valid.sum()), runs=len(np.unique(parent['run'][valid])),
        rejection_counts={str(k): int(v) for k, v in zip(*np.unique(reason, return_counts=True))},
        native_stored_beta_max_abs=float(np.max(abs(diff))) if len(diff) else None,
        native_stored_beta_exact=bool(np.all(diff == 0)),
        accepted_error_percentiles=np.percentile(error[valid], [0, 5, 50, 95, 100]) if valid.any() else [],
        reference_sha256=digest(REF / 'reference.npz'), protocol_sha256=digest(REF / 'protocol.json'),
        preparation_sha256=digest(REF / 'preparation.json'), output_freeze_sha256=p['output_freeze_sha256'],
        schema_amendment_sha256=digest(REF / 'schema_amendment.json'), parts=parts,
        new_RICH_values_read=True, no_output_based_selection=True, native_association_observed=True))
    print('PROSPECTIVE_RICH_REFERENCE', int(valid.sum()), 'of', n, 'runs', len(np.unique(parent['run'][valid])), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['amend', 'prepare', 'collect'])
    globals()[parser.parse_args().stage]()
