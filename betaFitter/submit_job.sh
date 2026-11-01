#!/bin/bash

Z=${1:-8}
fitOption=${2:-0}

arg4=${3:-1.0}
arg5=${4:-1.0}
arg6=${5:-all}
arg7=${6:-center}

condor_submit submit.sub -define "Z=${Z}" -define "fitOption=${fitOption}" \
    -define "arg4=${arg4}" -define "arg5=${arg5}" -define "arg6=${arg6}" \
    -define "arg7=${arg7}"
