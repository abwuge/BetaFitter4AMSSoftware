#!/usr/bin/env bash
set -e

source /home/ams/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh

if [[ $# -lt 2 ]]; then
    echo "Error: Need at least 2 parameters"
    echo "Usage: $0 inputFile outputFile [fitOption=0] [energyLossScale=1.0] [energyLossScaleMode=all]"
    exit 1
fi

inputFile="$1"
outputFile="$2"
fitOption="${3:-0}"
energyLossScale="${4:-1.0}"
energyLossScaleMode="${5:-all}"

# Get the script directory to locate the executable
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

"${SCRIPT_DIR}/build/bin/betaFitter" "$inputFile" "$outputFile" "$fitOption" "$energyLossScale" "$energyLossScaleMode"
