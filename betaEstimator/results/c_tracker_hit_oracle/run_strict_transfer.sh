#!/usr/bin/env bash
set -euo pipefail
cd /nas01/home/hxwu/AMSSoft/betaEstimator/results/c_tracker_hit_oracle
rm -rf tmva_strict_transfer ds_strictgeom_c12train_bin* ds_strictgeom_c13train_bin*
mkdir -p tmva_strict_transfer logs_strict_transfer
for spec in "strictgeom_c12train 0" "strictgeom_c13train 1"; do
  set -- $spec; stage=$1; label=$2; echo "START $stage"; pids=()
  for b in 0 1 2 3; do
    root -l -b -q "train_tracker_bdt_transfer.C(\"tracker_features.root\",\"tmva_strict_transfer\",\"${stage}\",\"features_strict_geometry.txt\",${label},${b},300)" > "logs_strict_transfer/${stage}_bin${b}.log" 2>&1 &
    pids+=($!)
  done
  fail=0; for pid in "${pids[@]}"; do if ! wait "$pid"; then fail=1; fi; done
  if [[ $fail -ne 0 ]]; then echo "FAILED $stage"; exit 1; fi
  echo "DONE $stage"
done
echo ALL_DONE
