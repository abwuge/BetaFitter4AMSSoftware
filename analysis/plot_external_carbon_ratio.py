#!/usr/bin/env python3
import json
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    out = ROOT / 'results/time_local_ratio/external_ratio'
    out.mkdir(exist_ok=True)
    reference = json.loads((out.parent / 'protocol.json').read_text())['reference']
    rows = list(zip(reference['x'], reference['ratio'], reference['error']))
    if not rows or not all(e > 0 and r > 0 and err >= 0 for e, r, err in rows):
        raise ValueError('Invalid external ratio points')
    with tempfile.TemporaryDirectory(prefix='external-carbon-ratio-', dir='/tmp') as temp:
        points = Path(temp) / 'points.txt'
        points.write_text(str(len(rows)) + '\n' + ''.join('%.17g %.17g %.17g\n' % row for row in rows))
        subprocess.run(['root', '-l', '-b', '-q', '%s("%s","%s")' %
            (ROOT / 'macro/plotExternalCarbonRatio.C', points, out)], check=True)
    subprocess.run(['pdfunite', str(out / 'external_ratio_energy.pdf'),
        str(out / 'external_ratio_energy_linear.pdf'), str(out / 'external_ratio.pdf')], check=True)


if __name__ == '__main__':
    main()
