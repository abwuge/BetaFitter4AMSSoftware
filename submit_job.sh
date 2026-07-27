#!/bin/bash

Z=${1:-8}
fitOption=${2:-0}
energyLossScale=${3:-1.0}
energyLossScaleMode=${4:-all}
referencePoint=${5:-center}

condor_submit submit.sub -define "Z=${Z}" -define "fitOption=${fitOption}" -define "energyLossScale=${energyLossScale}" -define "energyLossScaleMode=${energyLossScaleMode}" -define "referencePoint=${referencePoint}"
