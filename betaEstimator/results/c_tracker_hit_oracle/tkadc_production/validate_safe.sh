#!/usr/bin/env bash
set -o pipefail
source /home/ams/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh >/dev/null
prod=/nas01/home/hxwu/AMSSoft/betaEstimator/results/c_tracker_hit_oracle/tkadc_production
"$prod/analysis_amsd69n_ROOT5SLC6P8HPBNV11TKADCSAFE" "$1" "$2" 1 1000
status=$?
printf 'ANALYSIS_EXIT=%s\n' "$status"
if [[ $status -ne 3 ]]; then exit 20; fi
/usr/bin/python3 "$prod/validate_branches.py" "$2"
