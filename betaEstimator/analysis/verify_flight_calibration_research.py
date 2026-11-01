#!/usr/bin/env python3
"""Check research provenance and produce a compact result index."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/flight_calibration_research'
SCRIPTS={
    'input_inventory.json':'research_flight_calibration_inventory.py',
    'geometry_information.json':'research_flight_calibration_geometry.py',
    'linearized_closure.json':'research_flight_calibration_closure.py',
    'anchored_closure.json':'research_flight_calibration_anchor.py',
    'fixed_timing_closure.json':'research_flight_calibration_fixed_timing.py',
    'ends_response.json':'research_flight_calibration_ends.py'}


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    test=subprocess.run([sys.executable,str(ROOT/'tests/test_flight_calibration_core.py')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    (OUT/'unit_tests.log').write_text(test.stdout)
    if test.returncode:raise RuntimeError(test.stdout)
    core=digest(ROOT/'analysis/flight_calibration_core.py')
    results={}
    for name,script in SCRIPTS.items():
        d=json.loads((OUT/name).read_text())
        assert d['script_sha256']==digest(ROOT/'analysis'/script),(name,'script changed')
        if 'core_sha256' in d:assert d['core_sha256']==core
        if 'complete' in d:assert d['complete']
        results[name]=d
    for path in OUT.glob('*freeze.json'):
        d=json.loads(path.read_text())
        if 'core_sha256' in d:assert d['core_sha256']==core,path
    ends=results['ends_response.json'];freeze=json.loads((OUT/'ends_freeze.json').read_text())
    assert ends['freeze_sha256']==digest(OUT/'ends_freeze.json')
    assert not set(freeze['training_runs']).intersection(freeze['evaluation_runs'])
    assert ends['run_overlap']==0 and not ends['RICH_used'] and not ends['MC_used'] and not ends['beta_used']
    assert not ends['absolute_energy_calibration'] and not ends['beta_resolution_evaluation']
    assert ends['evaluated_common_hit_rows']==sum(r['pooled']['n'] for r in ends['rows'])
    for row in ends['rows']:assert row['pooled']['n']==row['time_conditioned']['n']
    fixed=results['fixed_timing_closure.json']
    fixed_rows=[]
    for r in fixed['rows']:
        a,b=r['metrics']['raw'],r['metrics']['energy_only_fit']
        fixed_rows.append(dict(domain=r['domain'],seed=r['seed'],w68_gain_pct=100*(1-b['w68']/a['w68']),
                               raw_bias=a['mean'],new_bias=b['mean'],rmse_improved=b['rmse']<a['rmse']))
    summary=dict(status='research_handoff_not_final_algorithm',
        completed=['nine-file measured-input inventory','100000-event flight-geometry information audit',
                   'linearized global-calibration primitive and six algebra tests',
                   'two-domain failed free-clock closure','global high-speed anchor repair test',
                   'eight fixed-clock energy-only closure trials','flight ADC-end ratio out-of-run pilot'],
        not_completed=['full raw-observable nonlinear likelihood','native Tracker-deposit supplement',
                       'flight-calibrated energy-loss model','AMS MC independent refit/closure',
                       'new blind RICH beta-resolution evaluation','global optimality or impossibility proof'],
        ends_response={k:ends[k] for k in ['measured_hit_rows','eligible_evaluation_hit_rows','evaluated_common_hit_rows','model_count','training_run_count','evaluation_run_count','run_overlap','rows']},
        fixed_timing_closure=fixed_rows,
        core_sha256=core,
        result_sha256={name:digest(OUT/name) for name in SCRIPTS},
        script_sha256={script:digest(ROOT/'analysis'/script) for script in SCRIPTS.values()})
    (OUT/'research_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    old_model_count=0
    for freeze_path in (ROOT/'results/time_local_beta/folds').glob('*/freeze.json'):
        old=json.loads(freeze_path.read_text())
        for block,expected in old['model_sha256'].items():
            assert digest(freeze_path.parent/('block_%s_model.json'%block))==expected
            old_model_count+=1
    assert old_model_count>0
    record=dict(passed=True,unit_test_count=6,script_and_core_hashes_match=True,
                ends_runs_disjoint=True,old_frozen_models_verified=old_model_count,
                scope='Provenance and implemented algebra only; not approval of the scientific model',
                verification_script_sha256=digest(__file__),summary_sha256=digest(OUT/'research_summary.json'))
    (OUT/'verification.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2),flush=True)
    print(json.dumps(fixed_rows,indent=2),flush=True)

if __name__=='__main__':main()
