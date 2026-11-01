#!/usr/bin/env python3
"""Verify the prospective flight experiment without rewriting any old freeze."""
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
import numpy as np
from research_flight_conditional import write
from research_conditional_mc import stamp
from prepare_prospective_flight import ROOT, OUT
from reference_prospective_flight_schema import REF, check as reference_check
from evaluate_prospective_flight import check as score_check
from bootstrap_prospective_flight import DEST, check as recalibration_check

BASE = ROOT / 'results/flight_calibration_research'
CONT = BASE / 'continuation_20261013'
SCORE = OUT / 'scoring'


def read(path):
    return json.loads(Path(path).read_text())


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def main():
    destination = OUT / 'verification.json'
    if destination.exists():
        raise FileExistsError('This prospective verification has already been saved')
    reference_check(); score_check(); recalibration_check()
    frozen = read(OUT / 'output_freeze.json'); audit = read(REF / 'audit.json')
    score = read(SCORE / 'evaluation.json'); recal = read(DEST / 'evaluation.json')
    endpoint = read(REF / 'endpoint_sensitivity/score_sensitivity.json')
    plane = read(REF / 'plane_audit.json'); native = read(OUT / 'native_verification.json')
    if not all(value['complete'] for value in [frozen, audit, score, recal, endpoint, plane, native]):
        raise ValueError('Incomplete prospective evidence')
    with np.load(OUT / 'measured.npz') as a, np.load(OUT / 'predictions_before_reference.npz') as pred, np.load(REF / 'reference.npz') as ref:
        for key in ['run', 'event', 'source', 'entry']:
            if not np.array_equal(a[key], pred[key]) or not np.array_equal(a[key], ref[key]):
                raise ValueError('Prospective global identity mismatch')
        mask = ref['valid']; ids = (a['run'].astype(np.uint64) << np.uint64(32)) | a['event'].astype(np.uint64)
        if len(np.unique(ids)) != len(ids) or int(mask.sum()) != score['reference_events'] or len(np.unique(a['run'][mask])) != score['runs']:
            raise ValueError('Prospective global cohort mismatch')
        with np.load(OUT / 'excluded_runs.npz') as old:
            if np.intersect1d(a['run'], old['run']).size:
                raise ValueError('Known historical run overlap')
        if any('rich' in key.lower() or 'truth' in key.lower() or key.startswith('mc_') for key in a.files):
            raise ValueError('Reference-like field found in measured input cache')
        measured_events = len(ids)
        output_coverage = {}
        for name in pred.files:
            if '_RICH_' in name and name.endswith('_u'):
                if not pred[name + '_valid'].all() or not np.isfinite(pred[name]).all() or np.any(pred[name] == 0):
                    raise ValueError('A frozen output does not cover every measured event')
                output_coverage[name] = len(pred[name])
        names = sorted(output_coverage)
        if len(names) != 90:
            raise ValueError('Unexpected prospective candidate count')
    for path, expected in [(OUT / 'predictions_before_reference.npz', frozen['prediction_sha256']),
        (REF / 'reference.npz', audit['reference_sha256']), (SCORE / 'bootstrap_metrics.npz', score['bootstrap_metrics_sha256']),
        (DEST / 'bootstrap_metrics.npz', recal['metrics_sha256'])]:
        if digest(path) != expected:
            raise ValueError('Prospective large output checksum mismatch: ' + str(path))
    if recal['replicates'] != 128 or recal['nested_combinations'] != 2048 or recal['coverage_failures'] or recal['numerical_training_failures']:
        raise ValueError('Recalibration replay incomplete or invalid')
    for name, expected in recal['parts_sha256'].items():
        path = DEST / 'parts' / name
        if digest(path) != expected or digest(path.with_suffix('.npz')) != read(path)['output_sha256']:
            raise ValueError('Prospective replay part changed')
    protected = read(CONT / 'inherited_status.json')['protected']
    for name, expected in protected.items():
        if digest(ROOT / name) != expected:
            raise ValueError('Inherited protected artifact changed: ' + name)
    old_count = 0; old_manifests = {}
    for path in sorted((ROOT / 'results/time_local_beta/folds').glob('*/freeze.json')):
        old_manifests[str(path.relative_to(ROOT))] = digest(path)
        for block, expected in read(path)['model_sha256'].items():
            if digest(path.parent / ('block_%s_model.json' % block)) != expected:
                raise ValueError('Old production model changed')
            old_count += 1
    if old_count != 999:
        raise ValueError('Unexpected old production-model count')
    training = ROOT / 'results/time_local_beta/training.npz'
    training_sha = digest(training)
    if training_sha != read(ROOT / 'results/time_local_beta/plan.json')['training_sha256']:
        raise ValueError('Original final training cache changed')
    prior = {'material_continuation_summary.json': '56d80e3bcd8f835d2adbb168a393e14902e1e6104ca5a2525a19d639fbd17e94',
        'material_continuation_verification.json': '5a3aa715425eda879b1af407d180c6337c569d3b05f8a1a7116ec8403994fc26'}
    for name, expected in prior.items():
        if digest(BASE / name) != expected:
            raise ValueError('Original handoff summary was rewritten')
    logs = OUT / 'verification_logs'; logs.mkdir(exist_ok=True)
    tests = {}
    for name in ['test_flight_calibration_core.py', 'test_flight_timing_redundancy.py', 'test_flight_mean_loss.py',
                 'test_constrained_joint_adc.py', 'test_continuous_rich_transport.py', 'test_prospective_evaluation.py']:
        result = subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', 'tests', '-p', name, '-v'],
            cwd=str(ROOT), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        log = logs / (name + '.log')
        if log.exists():
            raise FileExistsError('Do not overwrite an existing verification log')
        log.write_text(result.stdout)
        count = re.search(r'Ran (\d+) tests? in ', result.stdout)
        if result.returncode or count is None:
            raise RuntimeError('Prospective verification test failed: ' + result.stdout)
        tests[name] = dict(passed=True, tests=int(count.group(1)), log_sha256=digest(log))
    acceptance = []
    for index, (branch, prefix) in enumerate([('joint', 'newcalib_gaussianR_joint_fit_RICH_s1_u'),
                                           ('speed_gain_null', 'newcalib_gaussianR_speed_gain_null_fit_RICH_s1_u')]):
        intervals = recal['intervals']['beta']['joint'][branch + '_versus_ols']
        adjusted = intervals['two_primary_adjusted_interval']['w68_improvement']
        acceptance.append(dict(candidate=prefix, point=score['primary'][index]['result']['point'],
            evaluation_only_adjusted_interval=score['primary'][index]['result']['two_primary_adjusted_interval']['w68_improvement'],
            training_plus_evaluation_interval=intervals['interval']['w68_improvement'],
            training_plus_evaluation_adjusted_interval=adjusted,
            width_improvement_survives_training_and_two_primary_adjustment=bool(adjusted[0] > 0),
            calibration_interval_scope='Saved 256-node flight run refits and independent evaluation-run draws, not all upstream/systematic uncertainties'))
    queue = json.loads(subprocess.check_output(['condor_q', 'hxwu', '-json', '-attributes',
        'ClusterId,ProcId,JobStatus,RemoteHost']).decode() or '[]')
    record = dict(complete=True, created_utc=stamp(), numerical_provenance_and_integrity_passed=True,
        measured_events=measured_events, evaluated_events=score['reference_events'], evaluated_runs=score['runs'],
        predictions_frozen_before_reference=frozen['created_utc'] < read(REF / 'read_started.json')['created_utc'],
        output_candidates=len(output_coverage), recalibration_replicates=128, nested_evaluation_draws=2048,
        independent_evaluation_run_resamples=read(SCORE / 'protocol.json')['bootstrap']['replicates'],
        tests=tests, total_tests=sum(row['tests'] for row in tests.values()),
        protected_inherited_artifacts=len(protected), original_production_models_verified=old_count,
        original_production_freeze_manifests=old_manifests, original_final_training_sha256=training_sha,
        old_handoff_checksums=prior, primary_acceptance=acceptance,
        coordinate_equality_passed=plane['numerical_coordinate_match'], coordinate_sensitivity=endpoint['primary'],
        important_limitations=['Original nominal and native emission planes differ by up to 9.06 micrometers; equality check remains failed, with separately bounded endpoint sensitivity.',
            'The RICH reference cohort is new relative to enumerated persisted betaEstimator histories, not necessarily the upstream global TOF calibration.',
            'The two primary full-uncertainty beta-W68 tests do not pass positive lower bounds; no production replacement is accepted.',
            'Student response, two-group first-order loss, conditional selection and shared period-independent ADC calibration remain model assumptions.',
            'No global optimum or impossibility conclusion follows from this finite experiment.'],
        files_sha256={str(path.relative_to(ROOT)): digest(path) for path in [
            OUT / 'output_freeze.json', REF / 'audit.json', REF / 'schema_amendment.json', REF / 'plane_audit.json',
            REF / 'endpoint_sensitivity/evaluation.json', REF / 'endpoint_sensitivity/score_sensitivity.json',
            SCORE / 'protocol.json', SCORE / 'evaluation.json', DEST / 'protocol.json', DEST / 'point.json', DEST / 'evaluation.json',
            DEST / 'batch/python_runtime_amendment.json', DEST / 'batch/shared_runtime/manifest.json',
            DEST / 'batch/shared_runtime/probe_submission.json', DEST / 'batch/shared_runtime/remaining_submission.json']},
        active_user_scheduler_jobs=queue, source_sha256=digest(Path(__file__)), production_model_replaced=False)
    if not record['predictions_frozen_before_reference']:
        raise ValueError('Reference exposure predates prospective prediction freeze')
    write(destination, record)
    print('PROSPECTIVE_VERIFIED', json.dumps({key: record[key] for key in ['measured_events', 'evaluated_events', 'evaluated_runs',
        'total_tests', 'protected_inherited_artifacts', 'original_production_models_verified', 'primary_acceptance']}), flush=True)


if __name__ == '__main__':
    main()
