#!/usr/bin/env bash
set -e

source /home/ams/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh

if [[ $# -lt 2 ]]; then
    echo "Error: Need at least 2 parameters"
    echo "Usage: $0 inputFile outputFile [fitOption=0] [energyLossScale=1.0] [energyLossScaleMode=all|s1s2|s2] [referencePoint=center|before_tof] [globalZetaMin=0] [globalZetaMax=6]"
    echo "For fitOption=4, the fourth argument is betaMax instead of energyLossScale."
    exit 1
fi

inputFile="$1"
outputFile="$2"
fitOption="${3:-0}"
if [[ "$fitOption" == "4" ]]; then
    energyLossScale="${4:-0.9}"
else
    energyLossScale="${4:-1.0}"
fi
energyLossScaleMode="${5:-all}"
referencePoint="${6:-center}"
globalZetaMin="${7:-0}"
globalZetaMax="${8:-6}"

# Get the script directory to locate the executable
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

"${SCRIPT_DIR}/build/bin/betaFitter" "$inputFile" "$outputFile" "$fitOption" "$energyLossScale" "$energyLossScaleMode" "$referencePoint" "$globalZetaMin" "$globalZetaMax"
