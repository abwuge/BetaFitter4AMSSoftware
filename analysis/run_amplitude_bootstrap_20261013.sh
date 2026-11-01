#!/usr/bin/env bash
set -eo pipefail
cd /nas01/home/hxwu/AMSSoft/betaEstimator
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_NUM_THREADS=4
worker="$1"
base=results/flight_calibration_research/continuation_20261013/amplitude_noise/bootstrap/batch
pids=()
for offset in 0 32 64 96; do
  rep=$((worker+offset))
  python3 analysis/bootstrap_amplitude_noise_20261013.py worker --worker "$rep" --workers 128 > "$base/rep_${rep}.log" 2>&1 &
  pids+=("$!")
done
status=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then status=1; fi
done
printf 'AMPLITUDE_BOOT_WORKER_DONE %s status=%s host=%s\n' "$worker" "$status" "$(hostname)"
exit "$status"
