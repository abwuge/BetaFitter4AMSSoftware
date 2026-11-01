#!/usr/bin/env bash
set -euo pipefail
cd /nas01/home/hxwu/AMSSoft/betaEstimator/results/c_tracker_hit_oracle
rm -rf tmva_strict ds_strict_multifit_bin* ds_strict_fitdiag_bin* ds_strict_geometry_bin* ds_strict_morph_bin*
mkdir -p tmva_strict logs_strict
for spec in \
  "strict_multifit features_strict_multifit.txt" \
  "strict_fitdiag features_strict_fitdiag.txt" \
  "strict_geometry features_strict_geometry.txt" \
  "strict_morph features_strict_morph.txt"
do
  set -- $spec
  stage=$1; list=$2
  echo "START $stage"
  pids=()
  for b in 0 1 2 3; do
    root -l -b -q "train_tracker_bdt_list.C(\"tracker_features.root\",\"tmva_strict\",\"${stage}\",\"${list}\",${b},300)" > "logs_strict/${stage}_bin${b}.log" 2>&1 &
    pids+=($!)
  done
  fail=0
  for pid in "${pids[@]}"; do if ! wait "$pid"; then fail=1; fi; done
  if [[ $fail -ne 0 ]]; then
    echo "FAILED $stage"
    for b in 0 1 2 3; do tail -30 "logs_strict/${stage}_bin${b}.log" || true; done
    exit 1
  fi
  echo "DONE $stage"
done
echo ALL_DONE
