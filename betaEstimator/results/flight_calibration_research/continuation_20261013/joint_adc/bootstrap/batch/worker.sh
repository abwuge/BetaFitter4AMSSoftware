#!/usr/bin/env bash
set -eo pipefail
cd /nas01/home/hxwu/AMSSoft/betaEstimator
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
worker="$1"
pids=()
for offset in 0 32 64 96; do
  rep=$((worker+offset))
  python3 analysis/bootstrap_joint_adc_20261013.py worker --replicate "$rep" > "/nas01/home/hxwu/AMSSoft/betaEstimator/results/flight_calibration_research/continuation_20261013/joint_adc/bootstrap/batch/rep_${rep}.log" 2>&1 &
  pids+=("$!")
done
status=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then status=1; fi
done
printf 'JOINT_BOOTSTRAP_WORKER_DONE worker=%s status=%s host=%s\n' "$worker" "$status" "$(hostname)"
exit "$status"
