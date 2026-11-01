#!/usr/bin/env bash
set -euo pipefail
cd /nas01/home/hxwu/AMSSoft/betaEstimator/results/c_tracker_hit_oracle
rm -rf tmva_transfer ds_morph_c12train_bin* ds_morph_c13train_bin* ds_charge_c12train_bin* ds_charge_c13train_bin*
mkdir -p tmva_transfer logs_transfer
for spec in \
  "morph_c12train features_fitdiag_morph.txt 0" \
  "morph_c13train features_fitdiag_morph.txt 1" \
  "charge_c12train features_charge_probe.txt 0" \
  "charge_c13train features_charge_probe.txt 1"
do
  set -- $spec
  stage=$1; list=$2; label=$3
  echo "START $stage"
  pids=()
  for b in 0 1 2 3; do
    root -l -b -q "train_tracker_bdt_transfer.C(\"tracker_features.root\",\"tmva_transfer\",\"${stage}\",\"${list}\",${label},${b},300)" > "logs_transfer/${stage}_bin${b}.log" 2>&1 &
    pids+=($!)
  done
  fail=0
  for pid in "${pids[@]}"; do if ! wait "$pid"; then fail=1; fi; done
  if [[ $fail -ne 0 ]]; then
    echo "FAILED $stage"
    for b in 0 1 2 3; do tail -30 "logs_transfer/${stage}_bin${b}.log" || true; done
    exit 1
  fi
  echo "DONE $stage"
done
echo ALL_DONE
