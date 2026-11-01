#!/usr/bin/env bash
set -eo pipefail
cd /home/ams/hxwu/AMSSoft/betaEstimator
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
local_worker=${1:?worker index}
printf -v local_tasks 'results/time_local_beta/condor/fit_%03d.txt' "$local_worker"
python3 analysis/time_local_beta.py fit --batch "$local_tasks" --workers 32
