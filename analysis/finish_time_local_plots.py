#!/usr/bin/env python3
"""Finish plots and histogram checks after the main numerical pipeline."""
import json
from pathlib import Path
import subprocess
import sys
import time
import time_local_beta as t

while not (t.OUT/'verification.json').exists():
    log=Path('/tmp/betaEstimator-time-local-pipeline.log').read_text()
    if 'Traceback (most recent call last)' in log:
        raise RuntimeError('Main pipeline failed; inspect its log before plotting')
    time.sleep(30)
audit=json.loads((t.OUT/'verification.json').read_text())
if not audit['passed']:raise RuntimeError('Numerical verification failed')
for script in ('plot_time_local_beta.py','verify_time_local_beta.py'):
    subprocess.run([sys.executable,str(t.w.rt.ROOT/'analysis'/script)],check=True)
print('LOCAL_TIME_ALL_OUTPUTS_COMPLETE',flush=True)
