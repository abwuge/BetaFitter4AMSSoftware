#!/bin/bash
set -eo pipefail
cd /home/ams/hxwu/AMSSoft/betaEstimator
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 analysis/widen_material_mc.py worker --worker "$1"
