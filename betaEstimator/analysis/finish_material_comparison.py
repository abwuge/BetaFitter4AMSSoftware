#!/usr/bin/env python3
"""Export inverse-speed center comparisons and independent MC transport closure."""
import argparse
import json
import subprocess
from pathlib import Path
import numpy as np
import uproot
import rich_transport as rt
import evaluate_rich_transport as ev

OUT = rt.ROOT / 'results/material_comparison'


def load(path):
    with np.load(str(path), allow_pickle=False) as f:
        return {k: f[k] for k in f.files}


def stats(x):
    q = np.percentile(x, [16, 50, 84])
    return dict(n=len(x), mean=float(np.mean(x)), median=q[1], w68=(q[2]-q[0])/2)


def prepare(flight_path, mc_cache=None):
    OUT.mkdir(exist_ok=True)
    p = load(flight_path)
    reference = 1/p['rich_beta'] - p['delay_u_C12']
    new = p.get('ratio_corrected_u', p['local_corrected_u'])
    records = {}
    with uproot.recreate(str(OUT/'comparisons.root')) as f:
        for naf, name in [(True, 'NaF'), (False, 'AGL')]:
            s = p['reference_mask'] & (p['NaF'].astype(bool) == naf)
            s &= np.isfinite(reference) & np.isfinite(new)
            f['flight_'+name] = dict(R=p['R'][s], reference_beta=p['rich_beta'][s],
                u_reference=reference[s], u_raw=p['raw_u'][s], u_new=new[s],
                raw_residual=p['raw_u'][s]-reference[s], new_residual=new[s]-reference[s])
            records['flight_'+name] = dict(raw=stats(p['raw_u'][s]-reference[s]),
                new=stats(new[s]-reference[s]))
        if mc_cache is not None:
            cached = uproot.open(str(mc_cache/'comparisons.root'))
            old_metrics = json.loads((mc_cache/'metrics.json').read_text())
            for name in ['mc_NaF', 'mc_AGL']:
                f[name] = cached[name].arrays(library='np')
                records[name] = old_metrics['records'][name]
    if mc_cache is not None:
        save_metrics(flight_path, records)
        return
    with uproot.update(str(OUT/'comparisons.root')) as f:
        mc = load(rt.ROOT/'results/rich_transport/mc_predictions.npz')
        q = load(rt.ROOT/'results/rich_transport/mc_reference_v2_masks.npz')
        a = load(Path('/tmp/betaEstimator-rich-transport/mc_evaluation_input.npz'))
        for k in ['run', 'event']:
            if not np.array_equal(mc[k], a[k]):
                raise ValueError('MC association mismatch')
        geo = rt.MaterialMap().build_columns(a['tk_pos'], a['tk_dir'])
        nominal_u = np.sqrt(1+(11.17793/(6*mc['R']))**2)
        pred = rt.propagate(nominal_u, mc['R'], geo, mass_override=11.17793)
        delay = 1/pred['beta_planes'][:, -1] - nominal_u
        for naf, name in [(True, 'NaF'), (False, 'AGL')]:
            s = q['charge_consistent_reference'] & mc['truth_valid'] & pred['valid']
            s &= (mc['NaF'].astype(bool) == naf) & np.isfinite(mc['u_true'])
            reference = 1/mc['rich_beta'][s] - delay[s]
            residual = reference - mc['u_true'][s]
            f['mc_'+name] = dict(beta_mc=1/mc['u_true'][s], R=mc['R'][s],
                residual=residual, u_reference=reference, u_mc=mc['u_true'][s],
                isotope=mc['label'][s], delay_u=delay[s], source=mc['source'][s], run=mc['run'][s])
            records['mc_'+name] = stats(residual)
            records['mc_'+name]['reference_candidates'] = int((q['charge_consistent_reference'] &
                (mc['NaF'].astype(bool) == naf)).sum())
            records['mc_'+name]['isotopes'] = {str(int(label)): stats(residual[mc['label'][s] == label])
                for label in np.unique(mc['label'][s])}
    save_metrics(flight_path, records)


def save_metrics(flight_path, records):
    ev.dump(OUT/'metrics.json', dict(complete=True, flight_predictions=str(flight_path),
        flight_sha256=ev.digest(flight_path), definitions=dict(u='1/beta',
        u_material='fixed nominal C12 mean inverse-speed increment, L5 to RICH',
        u_reference='1/beta_RICH - u_material',
        closure='u_reference - u_MC(L5)',
        closure_scope='Combined material approximation, straggling and RICH measurement error'), records=records))
    print(json.dumps(records, default=lambda v: v.item() if isinstance(v,np.generic) else v), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--flight', type=Path, default=rt.ROOT/'results/time_local_beta/predictions.npz')
    parser.add_argument('--output', type=Path, default=OUT)
    parser.add_argument('--mc-cache', type=Path)
    args = parser.parse_args()
    OUT = args.output
    prepare(args.flight, args.mc_cache)
    subprocess.run(['root', '-l', '-b', '-q', '%s("%s",%s)' %
        (rt.ROOT/'macro/plotMaterialComparison.C', OUT, 'true' if 'ratio' in str(args.flight) else 'false')], check=True)
