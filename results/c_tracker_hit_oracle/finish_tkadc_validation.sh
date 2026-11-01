#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export ROOTSYS=/nas01/home/hxwu/.local
compiler=/cvmfs/sft.cern.ch/lcg/contrib/gcc/9.2.0/x86_64-centos7
export PATH="$ROOTSYS/bin:$compiler/bin:$PATH"
export LD_LIBRARY_PATH="$ROOTSYS/lib:$compiler/lib64:${LD_LIBRARY_PATH:-}"
cd /nas01/home/hxwu/AMSSoft/betaEstimator/results/c_tracker_hit_oracle
/usr/bin/python3 summarize_tkadc.py
"$ROOTSYS/bin/root" -l -b -q plot_tkadc_validation.C
echo TKADC_VALIDATION_DONE
