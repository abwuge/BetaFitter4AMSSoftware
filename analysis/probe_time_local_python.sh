#!/usr/bin/env bash
set -eo pipefail
hostname
python3 -c 'import sys,numpy,scipy,uproot; print(sys.version); print(numpy.__version__); print(scipy.__version__); print(uproot.__version__)'
cd /home/ams/hxwu/AMSSoft/betaEstimator
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -c 'import sys;sys.path.insert(0,"analysis");import time_local_beta; print("TRAINING_IMPORT_OK")'
