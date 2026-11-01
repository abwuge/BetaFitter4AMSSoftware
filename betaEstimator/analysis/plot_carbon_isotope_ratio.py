#!/usr/bin/env python3
"""Export frozen selected-sample carbon mixture ratios through interpreted ROOT."""
import csv
import io
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/carbon_isotope_ratio'
P = ROOT / 'results/flight_calibration_research/continuation_20261013/periodic_adc'
MASSES = (12 * .9314940954 - 6 * .0005109989461,
          13.00335483507 * .9314940954 - 6 * .0005109989461)
A = (12, 13)


def energy(r, i):
    return (math.hypot(6 * r, MASSES[i]) - MASSES[i]) / A[i]


def log_density_energy(t, e, i):
    kinetic = A[i] * e
    momentum = math.sqrt(kinetic * (kinetic + 2 * MASSES[i]))
    q = 6 / momentum
    z = math.log(q / .25)
    h = t[32] + t[33] * z
    logmix = -math.log1p(math.exp(h)) if i == 0 else -math.log1p(math.exp(-h))
    jacobian = 6 * A[i] * (kinetic + MASSES[i]) / momentum ** 3
    return t[30] * z + t[31] * z * z + logmix + math.log(jacobian)


def main():
    OUT.mkdir(exist_ok=True)
    selection = json.loads((P / 'conditional_noise/evaluation/measurement_comparison.json').read_text())['selected']['flight']
    specifications = [
        ('Selected NewCalib + ADC + noise', P / ('conditional_noise/flight/' + selection + '/freeze.json')),
        ('Native timing + Gaussian R', P / 'flight/native_gaussianR_step/freeze.json'),
        ('Native timing + Student-t R', P / 'flight/native_studentR5_step/freeze.json')]
    low = max(energy(2.5, i) for i in range(2))
    high = min(energy(6, i) for i in range(2))
    grid = [low * (high / low) ** (j / 200) for j in range(201)]
    models = []
    for label, source in specifications:
        d = json.loads(source.read_text())
        t = d['model']['theta']
        ratio = [math.exp(log_density_energy(t, e, 1) - log_density_energy(t, e, 0)) for e in grid]
        assert all(math.isfinite(v) and v > 0 for v in ratio)
        models.append(dict(label=label, source=str(source.relative_to(ROOT)),
            sha256=hashlib.sha256(source.read_bytes()).hexdigest(), theta30_34=t[30:35],
            energy_GeV_per_n=grid, ratio=ratio,
            checkpoints={str(e):math.exp(log_density_energy(t,e,1)-log_density_energy(t,e,0)) for e in [.65,.8,1.,1.5,1.9]}))
    old = json.loads((ROOT / 'results/flight_beta_wide/training.json').read_text())
    rows = []
    for block in old['blocks']:
        edges = block['edges']
        for j, f in enumerate(block['raw_spectrum']['fraction']):
            rlo, rhi = edges[j:j+2]
            center = math.sqrt(rlo*rhi)
            rows.append(dict(R_low_GV=rlo,R_high_GV=rhi,n=block['counts'][j],f13=f,
                ratio_same_R=f/(1-f),C12_equivalent_energy=energy(center,0),
                energy_low=energy(rlo,0),energy_high=energy(rhi,0),at_lower_bound=f<=.00100001))
    exported = (OUT/'reference_crdb.csv').read_text()
    records = list(csv.reader(io.StringIO(exported[exported.index('"13C/12C"'):])))
    measurements = []
    for record in records:
        assert record[0] == '13C/12C' and record[2] == 'EKN' and record[15] == '0'
        experiment = record[1]
        group = 2 if record[11] == '2025ApJ...993...81C' else (1 if experiment.startswith('Voyager') else 0)
        measurements.append(dict(experiment=experiment,group=group,
            energy=float(record[3]),low=float(record[4]),high=float(record[5]),ratio=float(record[6]),
            error_low=math.hypot(float(record[7]),float(record[9])),
            error_high=math.hypot(float(record[8]),float(record[10])),
            paper='https://ui.adsabs.harvard.edu/abs/'+record[11],
            distance_AU=float(record[13]),dates=record[14]))
    fraction_export = (OUT/'reference_fraction_crdb.csv').read_text()
    fraction_rows = list(csv.reader(io.StringIO(fraction_export[fraction_export.index('"13C/C"'):])))
    for r in fraction_rows:
        if r[11] != '1978ApJ...226..355B':
            continue
        f = float(r[6]); elo = math.hypot(float(r[7]),float(r[9])); ehi = math.hypot(float(r[8]),float(r[10]))
        ratio = f/(1-f)
        measurements.append(dict(experiment=r[1], group=3,energy=float(r[3]),low=float(r[4]),high=float(r[5]),ratio=ratio,
            error_low=ratio-(f-elo)/(1-f+elo),error_high=(f+ehi)/(1-f-ehi)-ratio,
            original_quantity='13C/C',original_fraction=f,original_error=[elo,ehi],
            paper='https://ui.adsabs.harvard.edu/abs/'+r[11],distance_AU=float(r[13]),dates=r[14]))
    references = dict(
        voyager=dict(ratio=.0629,error=.0033,energy_low=.048,energy_high=.126,
            source='https://doi.org/10.1086/176743',
            numerical_source='https://www.researchgate.net/publication/347818112_Carbon_Isotopes_in_the_PAMELA_Experiment',
            status='Voyager value quoted in Bogomolov et al. 2020'),
        pamela=dict(ratio=.06,energy_low=.1,energy_high=1.,
            source='https://doi.org/10.1088/1742-6596/1697/1/012004',
            status='Fixed analysis assumption in section 2; not an independent C13 abundance measurement'))
    references['CRDB'] = dict(fraction_query='https://lpsc.in2p3.fr/crdb/rest.php?num=13C&den=C&energy_type=EKN&format=csv', query='https://lpsc.in2p3.fr/crdb/rest.php?num=13C&den=12C&energy_type=EKN&format=csv', measurements=measurements, error='statistical and systematic errors added in quadrature', groups=['Near-Earth balloon/space', 'Voyager heliospheric', 'Voyager 2025 VLISM', 'Balloon 1977 (Buffington 1978)'])
    result = dict(complete=True,models=models,wide_training=rows,references=references,
        latest_definition='dN13/dT divided by dN12/dT at equal kinetic energy per nucleon; fitted selected-sample latent population',
        density_definition='dN_i/dq proportional to exp(theta30*z+theta31*z^2)*pi_i(q), z=ln(q/0.25), q=1/R',
        transform='q_i(T)=6/sqrt(A_i*T*(A_i*T+2*M_i)); dN_i/dT=(dN_i/dq)*abs(dq_i/dT)',
        measured_R_domain_GV=[2.5,6],common_energy_domain_GeV_per_n=[low,high],
        latest_supervision='Known C12/C13 masses; spectrum and isotope logodds fitted freely with detector response. No externally measured isotope ratio imposed.',
        wide_definition='Fixed same-R f13/(1-f13), with C12-equivalent energy axis; no equal-energy flux conversion',
        uncertainty='Frozen outputs supply no abundance confidence interval; no artificial binomial errors assigned',
        comparison_scope='Detector selection, relative isotope efficiency, fragmentation, rigidity migration and solar modulation remain in the comparison')
    (OUT/'ratios.json').write_text(json.dumps(result,indent=2)+'\n')
    with (OUT/'ratios.tsv').open('w') as f:
        f.write('energy_GeV_per_n\tselected_ratio\tnative_gaussianR_ratio\tnative_studentR_ratio\n')
        for j,e in enumerate(grid): f.write('\t'.join('%.12g'%v for v in [e]+[m['ratio'][j] for m in models])+'\n')
    with tempfile.TemporaryDirectory(prefix='carbon-ratio-',dir='/tmp') as tmp:
        path=Path(tmp)/'points.txt'
        with path.open('w') as f:
            f.write(str(len(grid))+'\n')
            for j,e in enumerate(grid): f.write(' '.join('%.17g'%v for v in [e]+[m['ratio'][j] for m in models])+'\n')
            f.write(str(len(rows))+'\n')
            for row in rows: f.write(' '.join(str(row[k]) for k in ['C12_equivalent_energy','energy_low','energy_high','ratio_same_R','n'])+' '+str(int(row['at_lower_bound']))+'\n')
            f.write(str(len(measurements))+'\n')
            for r in measurements:
                f.write(' '.join(str(r[k]) for k in ['group','energy','low','high','ratio','error_low','error_high'])+'\n')
        subprocess.run(['root','-l','-b','-q',str(ROOT/'macro/plotCarbonIsotopeRatio.C')+'("'+str(path)+'","'+str(OUT)+'")'],check=True)
    print(json.dumps(dict(output=str(OUT),energy_domain=[low,high],selected_checkpoints=models[0]['checkpoints']),indent=2))


if __name__ == '__main__':
    main()
