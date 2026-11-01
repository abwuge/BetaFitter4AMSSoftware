#!/usr/bin/env python3
"""Pin existing login NumPy/SciPy in shared storage without installing packages."""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path
import numpy
import scipy
from research_flight_conditional import write, digest
from research_conditional_mc import stamp
from bootstrap_prospective_flight import ROOT, DEST

BATCH = DEST / 'batch'
RUNTIME = BATCH / 'shared_runtime'


def queue(cluster):
    return json.loads(subprocess.check_output(['condor_q', str(cluster), '-json', '-attributes',
        'ClusterId,ProcId,JobStatus,Cmd,RemoteHost']).decode() or '[]')


def manifest(root):
    result = {}
    for path in sorted(root.rglob('*')):
        if path.is_file() and '__pycache__' not in path.parts and path.suffix not in ['.pyc', '.pyo']:
            result[str(path.relative_to(root))] = digest(path)
    return result


def prepare():
    if RUNTIME.exists():
        raise FileExistsError('Inspect an existing runtime snapshot rather than replacing it')
    if numpy.__version__ != '1.19.1' or scipy.__version__ != '0.18.1':
        raise ValueError('Unexpected login scientific runtime')
    RUNTIME.mkdir(); python = RUNTIME / 'python'; python.mkdir()
    sources = {'numpy': Path(numpy.__file__).parent, 'scipy': Path(scipy.__file__).parent,
               'numpy.libs': Path(numpy.__file__).parent.parent / 'numpy.libs'}
    versions = dict(numpy=numpy.__version__, scipy=scipy.__version__)
    records = {}
    for name, source in sources.items():
        before = manifest(source)
        shutil.copytree(str(source), str(python / name), ignore=shutil.ignore_patterns('__pycache__', '*.pyc', '*.pyo'))
        after = manifest(python / name)
        if before != after:
            raise ValueError('Scientific runtime snapshot differs: ' + name)
        records[name] = dict(source=str(source), files=before,
            total_bytes=sum(path.stat().st_size for path in (python / name).rglob('*') if path.is_file()))
    launcher = RUNTIME / 'worker.sh'
    launcher.write_text('''#!/usr/bin/env bash
set -eo pipefail
cd {root}
unset PYTHONHOME
export PATH=/usr/bin:/bin:$PATH
export PYTHONPATH=/home/ams/hxwu/.local/lib
/usr/bin/python3 -c 'import sys,numpy; print("NODE_DEFAULT_NUMPY",sys.executable,numpy.__version__,numpy.__file__,flush=True)'
export PYTHONPATH={python}:/home/ams/hxwu/.local/lib
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
/usr/bin/python3 -c 'import sys,numpy,scipy,awkward; assert numpy.__version__=="1.19.1" and scipy.__version__=="0.18.1"; assert hasattr(numpy,"AxisError") and awkward.__version__=="1.9.0"; print("PINNED_REPLAY_RUNTIME",sys.executable,numpy.__version__,numpy.__file__,scipy.__version__,awkward.__version__,flush=True)'
/usr/bin/python3 -m unittest discover -s tests -p test_prospective_evaluation.py -v
exec /usr/bin/python3 analysis/dispatch_prospective_recalibration.py worker --worker "$1"
'''.format(root=ROOT, python=python))
    launcher.chmod(0o755)
    env = dict(os.environ, PYTHONPATH=str(python) + ':/home/ams/hxwu/.local/lib',
        OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    env.pop('PYTHONHOME', None)
    output = subprocess.check_output(['/usr/bin/python3', '-c',
        'import numpy,scipy,awkward; print(numpy.__version__,numpy.__file__,scipy.__version__,scipy.__file__,awkward.__version__)'], env=env).decode()
    write(RUNTIME / 'manifest.json', dict(complete=True, created_utc=stamp(), versions=versions, packages=records,
        launcher_sha256=digest(launcher), source_sha256=digest(Path(__file__)), login_check=output,
        scientific_inference_protocol_sha256=digest(DEST / 'protocol.json'),
        method='Copy the already used login package trees and their bundled NumPy shared libraries to a new private shared-directory snapshot; prepend it only for this job. No package installation, system upgrade, monkey-patching, connector file export, or file download.',
        limitation='Operating-system libraries remain native; every worker runs the same import/version and four inference/metric tests before claiming work'))
    print('PINNED_REPLAY_RUNTIME_PREPARED', sum(row['total_bytes'] for row in records.values()), 'bytes', output.strip(), flush=True)


def apply(probe=False):
    stage = 'probe' if probe else 'remaining'
    destination = RUNTIME / (stage + '_submission.json')
    if destination.exists():
        raise FileExistsError('Pinned-runtime scheduler action already recorded')
    p = json.loads((RUNTIME / 'manifest.json').read_text())
    if not p['complete'] or digest(RUNTIME / 'worker.sh') != p['launcher_sha256']:
        raise ValueError('Pinned runtime launcher changed')
    cluster = json.loads((BATCH / 'submission.json').read_text())['cluster']
    before = queue(cluster)
    if probe:
        if any(row['JobStatus'] != 5 for row in before) or list((BATCH / 'claims').iterdir()):
            raise ValueError('A single-node runtime probe requires the untouched held queue')
        selected = [1]
    else:
        with (BATCH / 'worker_1.out').open() as stream:
            stdout = stream.read()
        if 'PINNED_REPLAY_RUNTIME' not in stdout or not list((BATCH / 'claims').glob('*/complete.json')):
            raise ValueError('First verify a real completed replay with the pinned runtime')
        selected = [row['ProcId'] for row in before if row['JobStatus'] == 5]
    selected_expression = ' || '.join('ProcId == %d' % number for number in selected)
    if not selected_expression:
        raise ValueError('No held workers to release')
    constraint = 'ClusterId == %d && JobStatus == 5 && (%s)' % (cluster, selected_expression)
    record = dict(created_utc=stamp(), stage=stage, cluster=cluster, worker_ids=selected, before=before,
        runtime_manifest_sha256=digest(RUNTIME / 'manifest.json'), old_runtime_amendment_sha256=digest(BATCH / 'python_runtime_amendment.json'),
        scientific_configuration_changed=False, previous_diagnosis_correction='Removing AMS setup alone did not resolve the issue: the worker operating systems do not provide the login machine NumPy version. The new launcher prints both the actual default and pinned versions.', complete=False)
    write(destination, record)
    record['qedit_response'] = subprocess.check_output(['condor_qedit', '-constraint', constraint, 'Cmd', json.dumps(str(RUNTIME / 'worker.sh'))], stderr=subprocess.STDOUT).decode()
    write(destination, record)
    record['release_response'] = subprocess.check_output(['condor_release', '-constraint', constraint], stderr=subprocess.STDOUT).decode()
    record['after'] = queue(cluster); record['complete'] = True
    write(destination, record); subprocess.check_call(['condor_reschedule'])
    print('PINNED_REPLAY_RUNTIME_RELEASED', stage, len(selected), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('stage', choices=['prepare', 'probe', 'remaining'])
    args = parser.parse_args()
    prepare() if args.stage == 'prepare' else apply(args.stage == 'probe')
