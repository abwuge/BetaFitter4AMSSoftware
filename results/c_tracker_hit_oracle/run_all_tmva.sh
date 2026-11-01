#!/usr/bin/env bash
set -euo pipefail
cd /nas01/home/hxwu/AMSSoft/betaEstimator/results/c_tracker_hit_oracle
rm -rf tmva_full ds_multifit_bin* ds_fitdiag_bin* ds_geometry_bin* ds_cluster_bin*
mkdir -p tmva_full logs
for spec in "multifit 194" "fitdiag 502" "geometry 611" "cluster 858"; do
  stage=${spec%% *}
  nvar=${spec##* }
  echo "START STAGE $stage nvar=$nvar"
  pids=()
  for b in 0 1 2 3; do
    root -l -b -q "train_tracker_bdt.C+(\"tracker_features.root\",\"tmva_full\",\"${stage}\",${nvar},${b},300)" > "logs/${stage}_bin${b}.log" 2>&1 &
    pids+=($!)
  done
  fail=0
  for pid in "${pids[@]}"; do
    if ! wait "$pid"; then fail=1; fi
  done
  if [[ $fail -ne 0 ]]; then
    echo "FAILED STAGE $stage"
    for b in 0 1 2 3; do tail -20 "logs/${stage}_bin${b}.log" || true; done
    exit 1
  fi
  echo "DONE STAGE $stage"
  for b in 0 1 2 3; do tail -3 "logs/${stage}_bin${b}.log"; done
done
echo ALL_DONE
