#!/usr/bin/env bash
set -eo pipefail
cd /nas01/home/hxwu/AMSSoft/betaEstimator
export PYTHONPATH=/home/ams/hxwu/.local/lib
unset PYTHONHOME
export PATH=/usr/bin:/bin:$PATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
/usr/bin/python3 -c 'import sys,numpy,awkward; assert numpy.__version__=="1.19.1" and hasattr(numpy,"AxisError"); assert awkward.__version__=="1.9.0"; print("REPLAY_RUNTIME",sys.executable,numpy.__version__,numpy.__file__,awkward.__version__,flush=True)'
exec /usr/bin/python3 analysis/dispatch_prospective_recalibration.py worker --worker "$1"
