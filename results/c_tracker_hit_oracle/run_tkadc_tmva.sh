#!/usr/bin/env bash
set -euo pipefail
cd /nas01/home/hxwu/AMSSoft/betaEstimator/results/c_tracker_hit_oracle
rm -rf tmva_tkadc ds_tkadc_base_bin* ds_tkadc_adcshape_bin* ds_tkadc_adcstrip_bin* ds_tkadc_morph_bin* ds_tkadc_shape_only_bin*
mkdir -p tmva_tkadc logs_tkadc
for spec in "tkadc_base features_tkadc_base.txt" "tkadc_adcshape features_tkadc_adcshape.txt" "tkadc_adcstrip features_tkadc_adcstrip.txt" "tkadc_morph features_tkadc_morph.txt" "tkadc_shape_only features_tkadc_shape_only.txt"
do
  set -- $spec
  stage=$1
  list=$2
  echo "START $stage"
  pids=()
  for b in 0 1 2 3; do
    root -l -b -q "train_tracker_bdt_list.C(\"tkadc_shape_features.root\",\"tmva_tkadc\",\"${stage}\",\"${list}\",${b},300)" > "logs_tkadc/${stage}_bin${b}.log" 2>&1 &
    pids+=($!)
  done
  fail=0
  for pid in "${pids[@]}"; do
    if ! wait "$pid"; then fail=1; fi
  done
  if [[ $fail -ne 0 ]]; then
    echo "FAILED $stage"
    for b in 0 1 2 3; do tail -40 "logs_tkadc/${stage}_bin${b}.log" || true; done
    exit 1
  fi
  echo "DONE $stage"
done
echo ALL_DONE
