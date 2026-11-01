#!/usr/bin/env bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export TKADC_WORKERS=1
cd /nas01/home/hxwu/AMSSoft/betaEstimator/results/c_tracker_hit_oracle
/usr/bin/python3 tkadc_production/audit_completed_production.py
/usr/bin/python3 extract_tkadc_shape_features.py
echo AUDIT_AND_EXTRACTION_DONE
