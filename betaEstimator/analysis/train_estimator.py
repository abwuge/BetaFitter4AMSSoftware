#!/usr/bin/env python3
"""Train on flight timing contrasts and freeze before any truth evaluation."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import uproot
from extract_inputs import CHARGES, load_arrays
from timing_model import clean, fit, score, reconstruct


def subset(d, mask):
    return {k: v[mask] for k, v in d.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cache', type=Path, required=True)
    ap.add_argument('--output', type=Path, default=Path('results/data_training'))
    ap.add_argument('--max-training-events', type=int, default=120000)
    args = ap.parse_args()
    if args.max_training_events < 500:
        ap.error('--max-training-events must be at least 500')
    args.output.mkdir(parents=True, exist_ok=True)
    models = args.output/'models'
    models.mkdir(exist_ok=True)
    audit = json.loads((args.cache/'extract_audit.json').read_text())
    summary = dict(training='measured TOF contrasts only', truth_access=False,
                   extraction=audit, max_training_events=args.max_training_events, charges={})
    summary['code_sha256'] = {name: hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
                             for name in ('train_estimator.py', 'timing_model.py', 'extract_inputs.py')}
    for z in CHARGES:
        tr = load_arrays(args.cache/('train_Z%d.npz' % z))
        va = load_arrays(args.cache/('validation_Z%d.npz' % z))
        if set(tr['run']) & set(va['run']):
            raise ValueError('Training/validation run overlap')
        tr = subset(tr, clean(tr))
        available = len(tr['run'])
        if available < 500:
            raise ValueError('Too few clean training events for Z=%d: %d' % (z, available))
        if available > args.max_training_events:
            index = np.random.RandomState(613+z).choice(available, args.max_training_events, replace=False)
            tr = subset(tr, np.sort(index))
        model = fit(tr, z)
        model_path = models/('model_Z%d.json' % z)
        model_path.write_text(json.dumps(model, indent=2))
        valid = clean(va)
        if valid.sum() < 100:
            raise ValueError('Too few held-out events for Z=%d' % z)
        heldout = score(subset(va, valid), model)
        u, applied = reconstruct(va, model)
        with np.errstate(invalid='ignore', divide='ignore'):
            beta = 1/u
        with uproot.recreate(str(args.output/('validation_Z%d.root' % z))) as output:
            output['validation'] = dict(run=va['run'], event=va['event'], sourceId=va['source_id'],
                Z=va['charge'], linearBeta=va['tof_betah'].astype(float),
                newBeta=beta, modelApplied=applied)
        summary['charges'][str(z)] = dict(available_training_entries=available,
            training=model['training'], validation=heldout, fallback_entries=int((~applied).sum()),
            model_sha256=hashlib.sha256(model_path.read_bytes()).hexdigest())
        print('Z=%d train=%d validation=%d heldout NLL improvement=%.5f' %
              (z, len(tr['run']), valid.sum(), heldout['nll_improvement']), flush=True)
    (args.output/'training_summary.json').write_text(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
