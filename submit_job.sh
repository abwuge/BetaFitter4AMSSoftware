#!/bin/bash

Z=${1:-8}
fitOption=${2:-0}

if [[ "$fitOption" == "4" ]]; then
    arg4=${3:-0.9}
    arg5=${4:-all}
    arg6=${5:-center}
    arg7=${6:-0}
    arg8=${7:-6}
    arg9=${8:-0}
    arg10=${9:-20}
else
    arg4=${3:-1.0}
    arg5=${4:-1.0}
    arg6=${5:-all}
    arg7=${6:-center}
    arg8=0
    arg9=0
    arg10=0
fi

condor_submit submit.sub -define "Z=${Z}" -define "fitOption=${fitOption}" \
    -define "arg4=${arg4}" -define "arg5=${arg5}" -define "arg6=${arg6}" \
    -define "arg7=${arg7}" -define "arg8=${arg8}" -define "arg9=${arg9}" \
    -define "arg10=${arg10}"
