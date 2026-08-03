#!/usr/bin/env bash
set -e

source /home/ams/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh

if [[ $# -lt 2 ]]; then
    echo "Error: Need at least 2 parameters"
    echo "Usage (reconstruction): $0 inputFile outputFile 0 [energyLossScale=1.0] [trackerEnergyLossScale=1.0] [energyLossScaleMode=all|s1s2|s2] [referencePoint=center|before_tof]"
    echo "Usage (joint fit): $0 inputFile outputFile 4 [betaMax=0.9] [energyLossScaleMode=all|s1s2|s2] [referencePoint=center|before_tof] [zetaMin=0] [zetaMax=6] [trackerScaleMin=0] [trackerScaleMax=20]"
    exit 1
fi

inputFile="$1"
outputFile="$2"
fitOption="${3:-0}"
if [[ "$fitOption" == "4" ]]; then
    betaMax="${4:-0.9}"
    energyLossScaleMode="${5:-all}"
    referencePoint="${6:-center}"
    globalZetaMin="${7:-0}"
    globalZetaMax="${8:-6}"
    trackerScaleMin="${9:-0}"
    trackerScaleMax="${10:-20}"
else
    energyLossScale="${4:-1.0}"
    trackerEnergyLossScale="${5:-1.0}"
    energyLossScaleMode="${6:-all}"
    referencePoint="${7:-center}"
fi

# Get the script directory to locate the executable
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [[ "$fitOption" == "4" ]]; then
    "${SCRIPT_DIR}/build/bin/betaFitter" "$inputFile" "$outputFile" "$fitOption" \
        "$betaMax" "$energyLossScaleMode" "$referencePoint" "$globalZetaMin" \
        "$globalZetaMax" "$trackerScaleMin" "$trackerScaleMax"
else
    "${SCRIPT_DIR}/build/bin/betaFitter" "$inputFile" "$outputFile" "$fitOption" \
        "$energyLossScale" "$trackerEnergyLossScale" "$energyLossScaleMode" "$referencePoint"
fi
