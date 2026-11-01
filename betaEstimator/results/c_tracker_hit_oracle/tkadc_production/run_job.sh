#!/usr/bin/env bash
set -o pipefail
source /home/ams/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh >/dev/null
input_root=$1
output_root=$2
nevents=${3:--1}
prod=/nas01/home/hxwu/AMSSoft/betaEstimator/results/c_tracker_hit_oracle/tkadc_production
binary="$prod/analysis_amsd69n_ROOT5SLC6P8HPBNV11TKADC"
tmp="${output_root}.part.${CONDOR_CLUSTER:-local}.${CONDOR_PROCESS:-0}.$$"
selected="${tmp}.selected"
trap 'rm -f "$tmp" "$selected"' EXIT
"$binary" "$input_root" "$tmp" 1 "$nevents"
status=$?
if [[ $status -ne 3 || ! -s "$tmp" ]]; then exit 20; fi
run=$(basename "$output_root" .root)
 /usr/bin/python3 "$prod/select_events.py" "$tmp" "$selected" "$prod/event_lists/$run.txt"
status=$?
if [[ $status -ne 0 || ! -s "$selected" ]]; then exit 21; fi
mv -f "$selected" "$output_root"
printf 'DONE %s\n' "$output_root"
