#!/usr/bin/env bash
set -eo pipefail
cd /nas01/home/hxwu/AMSSoft/betaEstimator
unset PYTHONHOME
export PATH=/usr/bin:/bin:$PATH
export PYTHONPATH=/nas01/home/hxwu/AMSSoft/betaEstimator/results/flight_calibration_research/continuation_20261013/prospective_flight/scoring/recalibration/batch/shared_runtime/python:/home/ams/hxwu/.local/lib
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
/usr/bin/python3 -c 'import numpy,scipy,awkward; assert numpy.__version__=="1.19.1" and scipy.__version__=="0.18.1"; assert awkward.__version__=="1.9.0"; print("PERIOD_RUNTIME",numpy.__version__,scipy.__version__,flush=True)'
/usr/bin/python3 -m unittest discover -s tests -p test_periodic_adc.py -v
exec /usr/bin/python3 analysis/dispatch_periodic_adc.py worker --index "$1"
