#!/usr/bin/env bash
set -eo pipefail
cd /home/ams/hxwu/AMSSoft/betaEstimator
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
ratio_worker=${1:?worker index}
printf -v ratio_tasks 'results/time_local_ratio/condor/fit_%03d.txt' "$ratio_worker"
python3 analysis/time_local_ratio.py fit --batch "$ratio_tasks" --workers 8
