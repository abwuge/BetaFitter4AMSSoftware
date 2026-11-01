#!/usr/bin/env bash
set -eo pipefail
cd /home/ams/hxwu/AMSSoft/betaEstimator
source /home/ams/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh
# Same shared GCC runtime used by the established TofCalib worker wrapper.
validation_gcc_runtime=/cvmfs/sft.cern.ch/lcg/contrib/gcc/5.3.0/x86_64-centos7-gcc53-opt/lib64
export LD_LIBRARY_PATH="$validation_gcc_runtime:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
validation_root="$PWD/results/tofcalib_validation"
worker_index=0
parallel_workers=32
single_source=""
runtime_probe=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        --worker) worker_index="$2"; shift 2 ;;
        --workers) shift 2 ;;
        --parallel) parallel_workers="$2"; shift 2 ;;
        --source) single_source="$2"; shift 2 ;;
        --runtime-probe) runtime_probe=true; shift ;;
        *) exit 2 ;;
    esac
done
native_binary="$validation_root/bin/supplement_light"
expected_hash=$(cat "$validation_root/bin/supplement_light.sha256")
actual_hash=$(sha256sum "$native_binary")
[[ "${actual_hash%% *}" == "$expected_hash" ]] || exit 3
if [[ "$runtime_probe" == true ]]; then
    printf -v source_name '%05d' "$single_source"
    task_directory="$validation_root/parts/$source_name"
    hostname
    "$native_binary" "$task_directory/evaluation_raw_files.txt" \
        "$task_directory/evaluation_targets.txt" "$validation_root/probe/el9_extra.root" \
        > "$validation_root/probe/el9_native.log" 2>&1
    exit $?
fi
run_source() {
    local task_directory="$validation_root/parts/$1"
    [[ -f "$task_directory/supplement.json" || -f "$task_directory/.native_complete" ]] && return 0
    "$native_binary" "$task_directory/evaluation_raw_files.txt" \
        "$task_directory/evaluation_targets.txt" "$task_directory/extra.partial.root" \
        > "$task_directory/evaluation_native.log" 2>&1 || return $?
    printf '%s\n' "$expected_hash" > "$task_directory/.native_complete.tmp"
    mv "$task_directory/.native_complete.tmp" "$task_directory/.native_complete"
    printf 'Native extraction complete: %s\n' "$1"
}
if [[ -n "$single_source" ]]; then
    printf -v source_name '%05d' "$single_source"
    run_source "$source_name"
else
    printf -v task_list '%s/condor/tasks_%03d.txt' "$validation_root" "$worker_index"
    run_partition() {
        local partition="$1" ordinal=0 source_name
        while IFS= read -r source_name; do
            if [[ -n "$source_name" ]] && (( ordinal % parallel_workers == partition )); then
                run_source "$source_name" || return $?
            fi
            ordinal=$((ordinal+1))
        done < "$task_list"
    }
    pids=()
    for ((partition=0; partition<parallel_workers; partition++)); do
        run_partition "$partition" &
        pids+=("$!")
    done
    failed=0
    for pid in "${pids[@]}"; do
        if wait "$pid"; then :; else failed=1; fi
    done
    exit "$failed"
fi
