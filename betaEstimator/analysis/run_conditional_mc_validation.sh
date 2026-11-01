#!/usr/bin/env bash
set -eo pipefail
cd /nas01/home/hxwu/AMSSoft/betaEstimator
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_NUM_THREADS=4
worker="$1"
python3 analysis/research_conditional_mc_reference.py run --worker "$worker" --workers 50 &
reference_pid=$!
set +e
python3 analysis/research_conditional_mc_bootstrap.py worker --worker "$worker" --workers 50
bootstrap_status=$?
wait "$reference_pid"
reference_status=$?
set -e
printf 'WORKER_DONE %s reference=%s bootstrap=%s host=%s\n' "$worker" "$reference_status" "$bootstrap_status" "$(hostname)"
if [ "$bootstrap_status" -ne 0 ] || [ "$reference_status" -ne 0 ]; then exit 1; fi
