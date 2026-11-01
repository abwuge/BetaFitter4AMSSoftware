#!/usr/bin/env python3
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import uproot

ROOT = Path(__file__).resolve().parents[1]


def main(wide_R=False):
    out = ROOT / 'results/time_local_ratio/material_closure'
    if wide_R:
        out = out.parent / 'material_closure_wide_R'
        if not json.loads((out / 'verification.json').read_text())['passed']:
            raise RuntimeError('Verify expanded MC evaluation first')
    verified = json.loads((out.parent / 'verification.json').read_text())
    prepared = json.loads((out.parent / 'material_closure/metrics.json').read_text())
    if not verified['passed'] or prepared['flight_sha256'] != verified['prediction_sha256']:
        raise RuntimeError('MC reference must match the verified evaluation snapshot')
    with uproot.open(str(out / 'comparisons.root')) as source:
        arrays = {rad: source['mc_' + rad].arrays(library='np') for rad in ['NaF', 'AGL']}
    with tempfile.TemporaryDirectory(prefix='material-mc-pairs-', dir='/tmp') as temp:
        points = Path(temp) / 'points.txt'
        with points.open('w') as stream:
            stream.write('2\n')
            for rad, a in arrays.items():
                x, y = a['beta_mc'], 1 / a['u_reference']
                if not (np.isfinite(x).all() and np.isfinite(y).all()):
                    raise ValueError('Non-finite MC reference pairs')
                stream.write('%s %d\n' % (rad, len(x)))
                for row in zip(x, y, a['isotope']):
                    stream.write('%.17g %.17g %d\n' % row)
        subprocess.run(['root', '-l', '-b', '-q', '%s("%s","%s",%.17g)' %
            (ROOT / 'macro/plotMaterialMCPairs.C', points, out, 200. if wide_R else 6.)], check=True)

    names = ['NaF_mc_beta_pairs', 'NaF_mc_residual_pairs', 'AGL_mc_beta_pairs',
        'AGL_mc_residual_pairs', 'AGL_mc_beta_pairs_zoom', 'AGL_mc_residual_pairs_zoom']
    subprocess.run(['pdfunite'] + [str(out / (name + '.pdf')) for name in names] +
        [str(out / 'material_mc_pairs.pdf')], check=True)


if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--wide-R',action='store_true')
    main(parser.parse_args().wide_R)
