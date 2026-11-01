#!/usr/bin/env bash
set -euo pipefail
cd /nas01/home/hxwu/AMSSoft/betaEstimator/results/c_tracker_hit_oracle
mkdir -p tmva_ablation logs_ablation
for spec in "otherhit features_otherhit.txt" "noother features_fitdiag_noother.txt" "morph features_fitdiag_morph.txt" "chargeprobe features_charge_probe.txt"; do
  stage=${spec%% *}
  list=${spec##* }
  echo "START $stage"
  pids=()
  for b in 0 1 2 3; do
    root -l -b -q "train_tracker_bdt_list.C+(\"tracker_features.root\",\"tmva_ablation\",\"${stage}\",\"${list}\",${b},300)" > "logs_ablation/${stage}_bin${b}.log" 2>&1 &
    pids+=($!)
  done
  fail=0
  for pid in "${pids[@]}"; do if ! wait "$pid"; then fail=1; fi; done
  if [[ $fail -ne 0 ]]; then
    echo "FAILED $stage"
    for b in 0 1 2 3; do tail -20 "logs_ablation/${stage}_bin${b}.log" || true; done
    exit 1
  fi
  echo "DONE $stage"
done
echo ALL_DONE
