#!/usr/bin/env python3
"""Final MC evaluation of frozen flight-trained models; truth is read only here."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import uproot
from extract_inputs import CHARGES, load_raw, parse_sources
from timing_model import reconstruct


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sources', type=Path, default=Path('config/sources.conf'))
    ap.add_argument('--models', type=Path, default=Path('results/data_training'))
    ap.add_argument('--beta-fitter', type=Path, default=Path('../betaFitter'))
    ap.add_argument('--output', type=Path, default=Path('results/data_trained_comparison'))
    args = ap.parse_args()
    frozen = json.loads((args.models/'training_summary.json').read_text())
    implementation = Path(__file__).parent/'timing_model.py'
    if hashlib.sha256(implementation.read_bytes()).hexdigest() != frozen['code_sha256']['timing_model.py']:
        raise ValueError('Timing implementation changed after training; retrain before evaluation')
    models, hashes = {}, {}
    for z in CHARGES:
        path = args.models/'models'/('model_Z%d.json' % z)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != frozen['charges'][str(z)]['model_sha256']:
            raise ValueError('Model changed after training: ' + str(path))
        hashes[str(z)] = digest
        models[z] = json.loads(path.read_text())
    parts = {z: [] for z in CHARGES}
    audit = []
    validation_index = 0
    for source_id, (role, source_z, path) in enumerate(parse_sources(args.sources), 1):
        if role != 'validation':
            continue
        validation_index += 1
        d, _ = load_raw(path)
        prediction = np.full(len(d['charge']), np.nan)
        applied = np.zeros(len(d['charge']), dtype=bool)
        for z in CHARGES:
            mask = d['charge'] == z
            u, ok = reconstruct({k: v[mask] for k,v in d.items()}, models[z])
            prediction[mask], applied[mask] = u, ok
        # Predictions are now fixed. The rest of this loop is evaluation only.
        with uproot.open(path) as f:
            truth = f['amstreea'].arrays(['mch','mmass','mevmom1'], library='np')
        truth = {k: v[d['source_entry']] for k,v in truth.items()}
        historical = (truth['mevmom1'][:,10] != -1000) & (truth['mmass'] > 0)
        part = args.beta_fitter/'results/hybrid_physical_20260913/validation/parts'/('part_%03d_Z%d.root' % (validation_index, source_z))
        with uproot.open(str(part)) as f:
            b = f['betaTree'].arrays(library='np')
        if len(b['Z']) != int(historical.sum()):
            raise ValueError('Historical source selection does not align: ' + str(part))
        if not np.array_equal(b['Z'], d['charge'][historical]):
            raise ValueError('Charge alignment failed: ' + str(part))
        if not np.array_equal(b['mcCharge'], truth['mch'][historical].astype(int)):
            raise ValueError('Truth charge alignment failed: ' + str(part))
        if not np.allclose(b['linearBeta'], d['tof_betah'][historical], rtol=0, atol=0, equal_nan=True):
            raise ValueError('Timing alignment failed: ' + str(part))
        p = truth['mevmom1'][historical,10].astype(float)
        mass = truth['mmass'][historical].astype(float)
        if not np.allclose(b['mcBeta'], p/np.sqrt(p*p+mass*mass), rtol=1e-6, atol=1e-7):
            raise ValueError('Truth beta alignment failed: ' + str(part))
        with np.errstate(invalid='ignore',divide='ignore'):
            new_beta = 1/prediction[historical]
        common = ((b['Z']==source_z) & (b['mcCharge']==source_z)
                  & np.isfinite(b['mcBeta']) & (b['mcBeta']>=.25) & (b['mcBeta']<=1))
        for key in ('linearBeta','nonlinearBeta','hybridBeta'):
            common &= np.isfinite(b[key]) & (b[key]>0)
        common &= np.isfinite(new_beta) & (new_beta>0)
        rows = {k: d[k][historical][common].astype(np.uint32) for k in ('run','event')}
        rows.update({k: b[k][common] for k in ('mcBeta','linearBeta','nonlinearBeta','hybridBeta','Z')})
        rows['sourceId'] = np.full(int(common.sum()),source_id,dtype=np.int32)
        rows['newBeta'] = new_beta[common]
        rows['modelApplied'] = applied[historical][common]
        parts[source_z].append(rows)
        audit.append(dict(path=path, source_id=source_id, aligned_rows=len(b['Z']),
                          common_entries=int(common.sum())))
        if validation_index % 16 == 0:
            print('evaluated %d validation sources' % validation_index, flush=True)
    args.output.mkdir(parents=True,exist_ok=True)
    summary = dict(model_hashes=hashes, training='flight data only', charges={}, alignment=audit)
    for z in CHARGES:
        if not parts[z]:
            raise ValueError('No evaluation sources for Z=%d' % z)
        rows = {k:np.concatenate([p[k] for p in parts[z]]) for k in parts[z][0]}
        dest = args.output/('comparison_Z%d.root' % z)
        with uproot.recreate(str(dest)) as f:
            f['comparison'] = rows
        n = len(rows['Z'])
        n06 = int(np.sum(rows['mcBeta']>=.6))
        summary['charges'][str(z)] = dict(common_entries=n,beta_ge_06_entries=n06,
            fallback_entries=int(np.sum(~rows['modelApplied'])),output=str(dest))
        print('Z=%d common=%d beta>=0.6=%d' % (z,n,n06),flush=True)
    (args.output/'comparison_summary.json').write_text(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
