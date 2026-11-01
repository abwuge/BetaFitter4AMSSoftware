#!/usr/bin/env python3
"""Check reference-plane coordinates without inspecting any beta residual."""
import json
from pathlib import Path
import numpy as np
from reference_prospective_flight_schema import ROOT, OUT, REF, check
from research_flight_conditional import load, write, digest
from research_conditional_mc import stamp


def main():
    check()
    destination = REF / 'plane_audit.json'
    if destination.exists():
        raise FileExistsError('Prospective reference-plane audit already exists')
    parent = load(REF / 'parent.npz')
    ref = load(REF / 'reference.npz')
    geometry = load(OUT / 'rich_geometry.npz')
    for key in ['run', 'event']:
        if not np.array_equal(parent[key], ref[key]) or not np.array_equal(geometry[key], ref[key]):
            raise ValueError('Reference-plane identity mismatch')
    mask = ref['valid'].astype(bool)
    actual = parent['rich_pos'][mask, 2].astype(float)
    target = geometry['target_z'][mask].astype(float)
    difference = actual - target
    if not np.isfinite(difference).all() or not len(difference):
        raise ValueError('Reference emission plane unavailable')
    match = bool(np.max(abs(difference)) < 1e-4)
    write(destination, dict(complete=True, created_utc=stamp(), events=int(mask.sum()),
        all_reference_NaF=bool(parent['rich_NaF'][mask].all()), actual_emission_z_cm=np.percentile(actual, [0, 5, 50, 95, 100]),
        propagation_target_z_cm=np.percentile(target, [0, 5, 50, 95, 100]),
        actual_minus_target_cm=np.percentile(difference, [0, 5, 50, 95, 100]),
        max_abs_difference_cm=float(np.max(abs(difference))), numerical_coordinate_match=match,
        numerical_tolerance_cm=1e-4, tolerance_purpose='Float-coordinate consistency, not a physical material-map accuracy claim',
        reference_sha256=digest(REF / 'reference.npz'), parent_sha256=digest(REF / 'parent.npz'),
        geometry_sha256=digest(OUT / 'rich_geometry.npz'), source_sha256=digest(Path(__file__)),
        reconstructed_beta_values_read=False, event_selection_changed=False))
    print('PROSPECTIVE_REFERENCE_PLANE', 'match', match, 'max_abs_cm', float(np.max(abs(difference))), flush=True)
    if not match:
        raise ValueError('Propagation endpoint differs from the observed reference plane; inspect before scoring')


if __name__ == '__main__':
    main()
