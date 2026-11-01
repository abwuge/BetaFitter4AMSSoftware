#!/usr/bin/env bash
set -o pipefail
source /home/ams/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh >/dev/null
input_root=$1
output_root=$2
nevents=${3:--1}
prod=/nas01/home/hxwu/AMSSoft/betaEstimator/results/c_tracker_hit_oracle/tkadc_production
binary="$prod/analysis_amsd69n_ROOT5SLC6P8HPBNV11TKADCSAFE"
tmp="${output_root}.part.${CONDOR_CLUSTER:-local}.${CONDOR_PROCESS:-0}.$$"
analysis_log="${tmp}.log"
trap 'rm -f "$tmp" "$analysis_log"' EXIT
"$binary" "$input_root" "$tmp" 1 "$nevents" 2>&1 | tee "$analysis_log"
status=${PIPESTATUS[0]}
if [[ $status -ne 3 || ! -s "$tmp" ]]; then exit 20; fi
run=$(basename "$output_root" .root)
/usr/bin/python3 "$prod/validate_full_output.py" "$tmp" "$analysis_log" "$prod/event_lists/$run.txt" "$nevents" "$input_root" || exit 21
mv -f "$tmp" "$output_root" || exit 22
printf 'DONE %s\n' "$output_root"
