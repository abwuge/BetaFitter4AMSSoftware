#!/usr/bin/env bash
set -eo pipefail
cd /home/ams/hxwu/AMSSoft/betaEstimator
source /home/ams/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh
local_gcc_runtime=/cvmfs/sft.cern.ch/lcg/contrib/gcc/5.3.0/x86_64-centos7-gcc53-opt/lib64
export LD_LIBRARY_PATH="$local_gcc_runtime:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
local_root="$PWD/results/time_local_beta"
local_binary="$PWD/results/tofcalib_validation/bin/supplement_light"
local_expected=$(cat "$PWD/results/tofcalib_validation/bin/supplement_light.sha256")
local_actual=$(sha256sum "$local_binary")
[[ "${local_actual%% *}" == "$local_expected" ]] || exit 3
local_worker=${1:?worker index}
local_parallel=${2:-32}
local_prefix=${3:-native}
printf -v local_list '%s/condor/%s_%03d.txt' "$local_root" "$local_prefix" "$local_worker"
run_partition() {
    local partition="$1" ordinal=0 source_name task_directory
    while IFS= read -r source_name; do
        if [[ -n "$source_name" ]] && (( ordinal % local_parallel == partition )); then
            task_directory="$local_root/parts/$source_name"
            if [[ -s "$task_directory/targets.txt" && ! -f "$task_directory/.native_complete" ]]; then
                "$local_binary" "$task_directory/raw_files.txt" "$task_directory/targets.txt" \
                    "$task_directory/extra.partial.root" > "$task_directory/native.log" 2>&1 || return $?
                mv "$task_directory/extra.partial.root" "$task_directory/extra.root"
                printf '%s\n' "$local_expected" > "$task_directory/.native_complete.tmp"
                mv "$task_directory/.native_complete.tmp" "$task_directory/.native_complete"
                printf 'Native training supplement complete: %s\n' "$source_name"
            fi
        fi
        ordinal=$((ordinal+1))
    done < "$local_list"
}
pids=()
for ((partition=0; partition<local_parallel; partition++)); do
    run_partition "$partition" &
    pids+=("$!")
done
failed=0
for pid in "${pids[@]}"; do
    if wait "$pid"; then :; else failed=1; fi
done
exit "$failed"
