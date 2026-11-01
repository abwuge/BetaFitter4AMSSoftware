#!/usr/bin/env bash
set -eo pipefail
cd "$(dirname "$0")/.."
source /nas01/home/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh
out=results/flight_calibration_research/continuation_20261013/prospective_flight
mkdir -p "$out/bin"
if [ -e "$out/bin/rich_reference" ]; then echo 'Reference binary exists; preserve it before rebuilding.' >&2; exit 2; fi
"$(root-config --cxx)" -O3 -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0 \
 -D_PGTRACK_ -DUSEGENFIT -D__ROOTSHAREDLIBRARY__ -D_IONL1PRESCALE_ -D_USELAYERRES_ \
 -I../amsd69n -isystem "$AMSWD/include" -isystem "$ROOTSYS/include" -isystem "$AMSWD/GenFit/include" \
 analysis/prospective_rich_reference.C -o "$out/bin/rich_reference" \
 -L"$AMSWD/lib/$(root-config --arch)5.34_7" -lntuple_slc6_PG \
 $(root-config --libs) -lMinuit -lTMVA -lNetx -lGeom -lMathMore -lEG -lTreePlayer -lMLP -lXMLIO \
 -L/cvmfs/ams.cern.ch/Offline/CERN/NagLib -lnaggcc64 \
 -L/cvmfs/ams.cern.ch/Offline/CERN/2005.gcc64.44/lib -lmathlib -lgfortran -lm -lgomp
sha256sum analysis/prospective_rich_reference.C "$out/bin/rich_reference"
