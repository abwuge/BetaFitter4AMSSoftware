#!/usr/bin/env python3
"""Remove an unnecessary AMS Python-path override from pure-Python replay jobs."""
import json
import subprocess
from pathlib import Path
from bootstrap_prospective_flight import ROOT, DEST
from research_flight_conditional import write, digest
from research_conditional_mc import stamp

BATCH = DEST / 'batch'


def queue(cluster):
    return json.loads(subprocess.check_output(['condor_q', str(cluster), '-json', '-attributes',
        'ClusterId,ProcId,JobStatus,Cmd,RemoteHost,HoldReason']).decode() or '[]')


def main():
    destination = BATCH / 'python_runtime_amendment.json'
    if destination.exists():
        raise FileExistsError('Pure-Python runtime amendment already recorded')
    cluster = json.loads((BATCH / 'submission.json').read_text())['cluster']
    before = queue(cluster)
    if any(row['JobStatus'] not in [1, 5] for row in before):
        raise ValueError('Inspect or finish active replay workers before changing an executable')
    if list((BATCH / 'claims').iterdir()) or list((DEST / 'parts').glob('*.json')):
        raise ValueError('This no-partial-result repair requires an untouched replay queue')
    launcher = BATCH / 'worker_python.sh'
    if launcher.exists():
        raise FileExistsError('Pure-Python launcher already exists')
    text = '''#!/usr/bin/env bash
set -eo pipefail
cd {root}
export PYTHONPATH=/home/ams/hxwu/.local/lib
unset PYTHONHOME
export PATH=/usr/bin:/bin:$PATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
/usr/bin/python3 -c 'import sys,numpy,awkward; assert numpy.__version__=="1.19.1" and hasattr(numpy,"AxisError"); assert awkward.__version__=="1.9.0"; print("REPLAY_RUNTIME",sys.executable,numpy.__version__,numpy.__file__,awkward.__version__,flush=True)'
exec /usr/bin/python3 analysis/dispatch_prospective_recalibration.py worker --worker "$1"
'''.format(root=ROOT)
    launcher.write_text(text); launcher.chmod(0o755)
    record = dict(created_utc=stamp(), cluster=cluster, before=before,
        reason='Sourcing the native AMS setup prepended an incompatible NumPy for Awkward import. These replay jobs do not need ROOT or AMS native libraries; use the already verified login interpreter and Python path instead.',
        original_launcher_sha256=digest(BATCH / 'worker.sh'), new_launcher_sha256=digest(launcher),
        source_sha256=digest(Path(__file__)), dispatch_protocol_sha256=digest(BATCH / 'protocol.json'),
        scientific_inference_protocol_sha256=digest(DEST / 'protocol.json'),
        calibration_and_reference_changed=False, software_packages_installed_or_upgraded=False,
        preexisting_replay_results=0, commands=[], complete=False)
    write(destination, record)
    if any(row['JobStatus'] == 1 for row in before):
        response = subprocess.check_output(['condor_hold', '-constraint', 'ClusterId == %d && JobStatus == 1' % cluster], stderr=subprocess.STDOUT).decode()
        record['commands'].append(dict(operation='hold unstarted workers', output=response)); write(destination, record)
    stopped = queue(cluster)
    if any(row['JobStatus'] != 5 for row in stopped):
        raise ValueError('Not all replay workers are held; inspect before changing the executable')
    response = subprocess.check_output(['condor_qedit', '-constraint', 'ClusterId == %d && JobStatus == 5' % cluster,
        'Cmd', json.dumps(str(launcher))], stderr=subprocess.STDOUT).decode()
    record['commands'].append(dict(operation='replace runtime launcher only', output=response)); write(destination, record)
    response = subprocess.check_output(['condor_release', '-constraint', 'ClusterId == %d && JobStatus == 5' % cluster], stderr=subprocess.STDOUT).decode()
    record['commands'].append(dict(operation='release same fifty workers', output=response))
    record['after'] = queue(cluster); record['complete'] = True
    write(destination, record); subprocess.check_call(['condor_reschedule'])
    print('PROSPECTIVE_REPLAY_RUNTIME_REPAIRED', len(before), 'unchanged worker identities', flush=True)


if __name__ == '__main__':
    main()
