#!/usr/bin/env python3
"""Show local-time width, bias, gain uncertainty and paired distributions."""
import json
from pathlib import Path
import numpy as np
import time_local_beta as t


def plot_distributions(out, r, p):
    """Bridge the existing NumPy prediction cache to the ROOT plotting macro."""
    import subprocess
    import tempfile
    distribution = out / 'distributions'
    distribution.mkdir(exist_ok=True)
    samples_all = [p[k] - p['rich_beta'] for k in ('raw_beta_rich', 'local_beta_rich')]
    histograms = {}
    pages = []
    for axis, records in r['rows'].items():
        values = p['R'] if axis == 'rigidity' else p['rich_beta']
        axis_records = []
        for index, row in enumerate(records):
            mask = p['reference_mask'] & (values >= row['low']) & (values < row['high'])
            samples = [v[mask] for v in samples_all]
            n = int(mask.sum())
            note = dict(low=row['low'], high=row['high'], n=n, views={})
            fit_inputs = []
            for sample in samples:
                q25, median, q75 = np.percentile(sample, [25, 50, 75])
                sigma = max((q75-q25)*.741301109252801, 1e-7)
                flo, fhi = median-3*sigma, median+3*sigma
                counts, _ = np.histogram(sample, np.linspace(flo, fhi, 121))
                fit_inputs.append((flo, fhi, median, sigma, counts))
            for view in ('core', 'full'):
                merged = np.concatenate(samples)
                low, high = (merged.min(), merged.max()) if view == 'full' else np.percentile(merged, [.5, 99.5])
                if high <= low: low, high = low-.001, high+.001
                edges = np.linspace(low, high, 121)
                counts = [np.histogram(sample, edges)[0] for sample in samples]
                outside = [int(n-h.sum()) for h in counts]
                note['views'][view] = dict(edges=edges, counts=counts, outside=outside)
                pages.append((axis, view, index, row, n, low, high, counts, fit_inputs))
            axis_records.append(note)
        histograms[axis] = axis_records
    pages.sort(key=lambda item: (item[0], item[1], item[2]))
    with tempfile.TemporaryDirectory(prefix='beta-distributions-', dir='/tmp') as scratch:
        source = Path(scratch) / 'histograms.txt'
        with source.open('w') as stream:
            stream.write(str(len(pages))+'\n')
            for axis, view, index, row, n, low, high, counts, fit_inputs in pages:
                stream.write('%s %s %d %.17g %.17g %d 120 %.17g %.17g\n' % (axis, view, index, row['low'], row['high'], n, low, high))
                for h in counts: stream.write(' '.join(map(str, h))+'\n')
                for flo, fhi, median, sigma, h in fit_inputs:
                    stream.write('120 %.17g %.17g %.17g %.17g\n' % (flo, fhi, median, sigma))
                    stream.write(' '.join(map(str, h))+'\n')
        macro = Path(__file__).resolve().parents[1] / 'macro/plotTimeLocalDistributions.C'
        subprocess.run(['root', '-l', '-b', '-q', '%s("%s","%s")' % (macro, source, distribution)], check=True)
    t.w.ev.dump(distribution/'histograms.json', dict(complete=True, axes=histograms, methods=['raw', 'local']))


def plot_trends(out, r):
    import subprocess
    import tempfile
    from datetime import datetime, timezone
    dest = out/'figures'
    dest.mkdir(exist_ok=True)
    groups = []
    for axis, records in r['rows'].items():
        rows = [row for row in records if row['n'] >= 20]
        bounds = (2.5, 200) if axis == 'rigidity' else (records[0]['low'], records[-1]['high'])
        groups.append((axis, rows, bounds))
        if axis == 'beta': groups.append(('beta_zoom', rows, (.975, 1.004)))
    for axis, records in [('time', r['period_rows']), ('yearly', r['year_rows'])]:
        rows = []
        for record in records:
            if 'local' not in record: continue
            row = dict(record['local'])
            row['x'] = (record['start']+record['stop'])/2 if axis == 'time' else datetime(record['year'], 7, 1, tzinfo=timezone.utc).timestamp()
            row['bootstrap'] = record.get('bootstrap', {'beta_gain_pct_q025_q50_q975': [row['gain_beta_pct']]*3})
            rows.append(row)
        groups.append((axis, rows, (rows[0]['x']-180*86400, rows[-1]['x']+180*86400)))
    with tempfile.TemporaryDirectory(prefix='beta-trends-', dir='/tmp') as scratch:
        source = Path(scratch)/'points.txt'
        with source.open('w') as stream:
            stream.write(str(len(groups))+'\n')
            for axis, rows, bounds in groups:
                stream.write('%s %d %.17g %.17g\n' % (axis, len(rows), *bounds))
                for row in rows:
                    interval = row['bootstrap']['beta_gain_pct_q025_q50_q975']
                    values = [row['x']] + [row[key]['w68']*1000 for key in ('baseline_beta','new_beta')] + [row[key]['mean']*1000 for key in ('baseline_beta','new_beta')] + [row['gain_beta_pct'], interval[0], interval[2]]
                    assert np.isfinite(values).all()
                    stream.write(' '.join('%.17g' % value for value in values)+'\n')
        macro = Path(__file__).resolve().parents[1]/'macro/plotTimeLocalTrends.C'
        subprocess.run(['root', '-l', '-b', '-q', '%s("%s","%s")' % (macro, source, dest)], check=True)
    for axis in ('time', 'yearly'):
        (dest/(axis+'_comparison.pdf')).unlink()
    for extension in ('pdf', 'png'):
        (dest/('beta_comparison_zoom.'+extension)).write_bytes((dest/('beta_zoom_comparison.'+extension)).read_bytes())
        (dest/('beta_zoom_comparison.'+extension)).unlink()


def main():
    out = t.OUT
    r = json.loads((out/'results.json').read_text())
    plot_trends(out, r)
    p = t.load(out/'predictions.npz')
    plot_distributions(out, r, p)
    print('LOCAL_TIME_PLOTS_COMPLETE',flush=True)


if __name__=='__main__':
    import sys
    if '--trends-only' in sys.argv:
        plot_trends(t.OUT, json.loads((t.OUT/'results.json').read_text()))
    elif '--distributions-only' in sys.argv:
        plot_distributions(t.OUT, json.loads((t.OUT/'results.json').read_text()), t.load(t.OUT/'predictions.npz'))
    else:
        main()
