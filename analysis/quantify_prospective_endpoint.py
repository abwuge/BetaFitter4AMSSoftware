#!/usr/bin/env python3
"""Quantify the observed sub-millimeter reference-endpoint discrepancy.

Keep every original frozen prediction unchanged. Re-integrate only the last
material cell to each observed emission z as a separately labelled sensitivity.
No beta residual is evaluated or used to select an endpoint or candidate.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import rich_transport as material
import continuous_rich_transport as forward
from reference_prospective_flight_schema import ROOT, OUT, REF, check
from research_flight_conditional import load, write, digest
from research_conditional_mc import stamp

DEST = REF / 'endpoint_sensitivity'


def plan():
    check()
    DEST.mkdir(exist_ok=True)
    if (DEST / 'protocol.json').exists():
        raise FileExistsError('Endpoint sensitivity protocol already recorded')
    audit = json.loads((REF / 'plane_audit.json').read_text())
    if audit['numerical_coordinate_match']:
        raise ValueError('No failed coordinate audit to investigate')
    if (OUT / 'scoring/evaluation.json').exists():
        raise ValueError('Record endpoint sensitivity before scoring beta residuals')
    write(DEST / 'protocol.json', dict(created_utc=stamp(), source_sha256=digest(Path(__file__)),
        original_audit_sha256=digest(REF / 'plane_audit.json'),
        prediction_freeze_sha256=digest(OUT / 'output_freeze.json'),
        reference_sha256=digest(REF / 'reference.npz'), parent_sha256=digest(REF / 'parent.npz'),
        geometry_sha256=digest(OUT / 'rich_geometry.npz'), map_sha256=digest(material.MAP),
        material_implementation_sha256=digest(Path(material.__file__)),
        forward_implementation_sha256=digest(Path(forward.__file__)),
        reason='Observed native RICH emission z differs from the nominal -75.552 cm endpoint by up to about 9 micrometers; the original 1e-4 cm numerical-coordinate check failed and is not retrospectively changed.',
        scope='Same frozen reference cohort and all ninety original candidate/material variants; no refit or output-based event removal',
        method='Recompute the final elemental material column to each stored emission z with the identical five-point sampling and float32 storage, then replay continuous propagation from the already frozen L5 inverse speed.',
        outputs='Original predictions remain the nominal primary outputs. Observed-z predictions are separately saved as a coordinate sensitivity, never chosen using reference residuals.',
        deterministic_width_bound='For a fixed cohort with max absolute beta-output change epsilon, every empirical beta-residual quantile changes by at most epsilon, and W68=(Q84-Q16)/2 changes by at most epsilon. This bound does not require reading reference beta.',
        resolution_claim='This measures endpoint-choice sensitivity within the existing material/propagation model, not the accuracy of the material map, RICH alignment, or energy-loss law.',
        reconstructed_beta_residuals_evaluated=False, selection_or_primary_definition_changed=False))
    print('PROSPECTIVE_ENDPOINT_SENSITIVITY_PROTOCOL', flush=True)


def run():
    check()
    p = json.loads((DEST / 'protocol.json').read_text())
    checks = [(Path(__file__), 'source_sha256'), (REF / 'plane_audit.json', 'original_audit_sha256'),
        (OUT / 'output_freeze.json', 'prediction_freeze_sha256'), (REF / 'reference.npz', 'reference_sha256'),
        (REF / 'parent.npz', 'parent_sha256'), (OUT / 'rich_geometry.npz', 'geometry_sha256'),
        (material.MAP, 'map_sha256'), (Path(material.__file__), 'material_implementation_sha256'),
        (Path(forward.__file__), 'forward_implementation_sha256')]
    for path, key in checks:
        if digest(path) != p[key]:
            raise ValueError('Endpoint sensitivity input changed: ' + key)
    if (DEST / 'evaluation.json').exists() or (DEST / 'predictions.npz').exists():
        raise FileExistsError('Endpoint sensitivity already saved')
    with np.load(REF / 'reference.npz') as f:
        mask = f['valid'].copy(); identities = {key: f[key][mask] for key in ['run', 'event']}
    parent = load(REF / 'parent.npz'); measured = load(OUT / 'measured.npz')
    old_geometry = load(OUT / 'rich_geometry.npz')
    geometry = forward.subset_geometry(old_geometry, mask)
    actual = parent['rich_pos'][mask, 2].astype(float)
    high = float(geometry['faces'][-2])
    if np.any(actual >= high) or not np.all(geometry['NaF_predicted']):
        raise ValueError('Observed emission endpoint is outside the last NaF cell')
    mapper = material.MaterialMap()
    new_column, inside = mapper.column_density(measured['tk_pos'][mask], measured['tk_dir'][mask],
        np.full(len(actual), high), actual, samples=5)
    if not inside.all() or not np.isfinite(new_column).all() or np.any(new_column < 0):
        raise ValueError('Observed endpoint has invalid material column')
    previous = geometry['column'][:, -1].copy()
    geometry['column'] = geometry['column'].copy()
    geometry['column'][:, -1] = new_column.astype(np.float32)
    geometry['target_z'] = actual
    np.savez_compressed(DEST / 'geometry.npz', **geometry)
    predicted = load(OUT / 'predictions_before_reference.npz')
    names = sorted(key for key in predicted if '_RICH_' in key and key.endswith('_u'))
    values = dict(identities); rows = {}
    for name in names:
        center_name, ending = name.split('_RICH_')
        scale_name = ending[:-2]
        scale = {'s1': 1., 's0p5': .5, 's2': 2.}[scale_name]
        u = predicted[center_name + '_L5_u'][mask]
        result = forward.propagate(u, measured['R'][mask], geometry, material_scale=scale)
        if not result['valid'].all() or not np.isfinite(result['u']).all() or np.any(result['u'] == 0):
            raise ValueError('An observed-endpoint output failed; do not filter the cohort')
        values[name] = result['u']
        old = predicted[name][mask]
        du = result['u'] - old
        db = 1 / result['u'] - 1 / old
        rows[name] = dict(max_abs_du=float(np.max(abs(du))), max_abs_dbeta=float(np.max(abs(db))),
            RMS_dbeta=float(np.sqrt(np.mean(db ** 2))), beta_width_absolute_bound=float(np.max(abs(db))))
    np.savez_compressed(DEST / 'predictions.npz', **values)
    write(DEST / 'evaluation.json', dict(complete=True, created_utc=stamp(), events=len(actual), candidates=len(names),
        original_coordinate_match=False, original_numerical_tolerance_changed=False,
        last_cell_high_z_cm=high, signed_last_cell_column_change_gcm2=np.percentile(
            geometry['column'][:, -1].sum(1) - previous.sum(1), [0, 5, 50, 95, 100]),
        max_abs_du=max(row['max_abs_du'] for row in rows.values()),
        max_abs_dbeta=max(row['max_abs_dbeta'] for row in rows.values()), branches=rows,
        original_freeze_sha256=digest(OUT / 'output_freeze.json'), reference_sha256=p['reference_sha256'],
        protocol_sha256=digest(DEST / 'protocol.json'), geometry_sha256=digest(DEST / 'geometry.npz'),
        predictions_sha256=digest(DEST / 'predictions.npz'), primary_predictions_changed=False,
        reconstructed_beta_residuals_evaluated=False, no_event_removal=True))
    print('PROSPECTIVE_ENDPOINT_SENSITIVITY', len(actual), len(names),
        'max_abs_dbeta', max(row['max_abs_dbeta'] for row in rows.values()), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('stage', choices=['plan', 'run'])
    globals()[parser.parse_args().stage]()
