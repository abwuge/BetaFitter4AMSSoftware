#!/usr/bin/env bash
set -eo pipefail
cd /nas01/home/hxwu/AMSSoft/betaEstimator
source /nas01/home/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec python3 analysis/dispatch_prospective_reference.py worker --worker "$1"
