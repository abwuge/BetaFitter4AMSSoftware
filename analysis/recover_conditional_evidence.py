#!/usr/bin/env python3
"""Recheck the saved results behind the interrupted user-facing deliveries."""
import json
from pathlib import Path
from research_flight_conditional import ROOT, digest, write
from research_conditional_mc import stamp

BASE = ROOT / 'results/flight_calibration_research'
CONT = BASE / 'continuation_20261013'
DEST = CONT / 'periodic_adc/recovered_evidence.json'


def main():
    if DEST.exists():
        raise FileExistsError('Recovered evidence was already recorded')
    original = json.loads((CONT / 'inherited_status.json').read_text())
    verified = json.loads((CONT / 'prospective_flight/verification.json').read_text())
    checks = dict(original['protected']); checks.update(verified['files_sha256'])
    for name, expected in checks.items():
        if digest(ROOT / name) != expected:
            raise ValueError('Saved evidence changed: ' + name)
    result = dict(complete=True, created_utc=stamp(), inherited_artifacts_rechecked=len(checks),
        source_sha256=digest(Path(__file__)), old_delivery_status='Results are present on the server; the visible transcript contained failed assistant turns. The failure mechanism is unconfirmed.',
        source_results={}, flight={}, additional_mc={}, prior_constrained_training={})
    flight = json.loads((CONT / 'prospective_flight/scoring/evaluation.json').read_text())
    result['flight'] = dict(measured_events=flight['measured_events'], reference_events=flight['reference_events'], runs=flight['runs'],
        primary_full_uncertainty=verified['primary_acceptance'], limitations=verified['important_limitations'], beta_metrics={})
    for branch in ['ols','time','joint','speed_null','speed_gain_null']:
        name = 'newcalib_gaussianR_' + branch + '_fit_RICH_s1_u'
        result['flight']['beta_metrics'][branch] = flight['metrics']['beta'][name]
        print('RECOVERED_FLIGHT_BETA', branch, json.dumps(flight['metrics']['beta'][name]), flush=True)
    mc = json.loads((CONT / 'additional_mc/scoring/evaluation.json').read_text())
    result['additional_mc'] = dict(events=mc['events'], runs=mc['runs'], beta={})
    for name, row in mc['table']['beta']['branches'].items():
        if name == 'native_u' or '_joint_' in name or name.endswith('time_transport_u'):
            result['additional_mc']['beta'][name] = row
            print('RECOVERED_MC_BETA', name, 'W68', row['metrics']['w68'], 'gain', row['versus_native']['w68_improvement'],
                'adjusted_interval', row['candidate_adjusted_interval']['w68_improvement'], flush=True)
    con = json.loads((CONT / 'constrained_adc/retraining_evaluation.json').read_text())
    result['prior_constrained_training'] = {key: con[key] for key in ['events','runs','replicates','nested_combinations','numerically_accepted','convergence_failures','intervals']}
    for relative in ['inherited_status.json','prospective_flight/verification.json','prospective_flight/scoring/evaluation.json',
                     'additional_mc/scoring/evaluation.json','constrained_adc/retraining_evaluation.json']:
        result['source_results'][relative] = digest(CONT / relative)
    result['scope'] = 'Rechecked inherited evidence, not newly performed calibration or scientific acceptance'
    write(DEST, result)
    print('RECOVERED_EVIDENCE_VERIFIED', len(checks), 'artifacts', digest(DEST), flush=True)


if __name__ == '__main__':
    main()
