#!/usr/bin/env bash
set -eo pipefail
cd "$(dirname "$0")/.."
source /nas01/home/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh
"$(root-config --cxx)" -O3 -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0 \
 -D_PGTRACK_ -DUSEGENFIT -D__ROOTSHAREDLIBRARY__ -D_IONL1PRESCALE_ -D_USELAYERRES_ \
 -I../amsd69n -isystem "$AMSWD/include" -isystem "$ROOTSYS/include" -isystem "$AMSWD/GenFit/include" \
 analysis/conditional_mc_reference_parent.C -o results/flight_calibration_research/conditional_mc/reference/extract_parent \
 -L"$AMSWD/lib/$(root-config --arch)5.34_7" -lntuple_slc6_PG \
 $(root-config --libs) -lMinuit -lTMVA -lNetx -lGeom -lMathMore -lEG -lTreePlayer -lMLP -lXMLIO \
 -L/cvmfs/ams.cern.ch/Offline/CERN/NagLib -lnaggcc64 \
 -L/cvmfs/ams.cern.ch/Offline/CERN/2005.gcc64.44/lib -lmathlib -lgfortran -lm -lgomp
sha256sum analysis/conditional_mc_reference_parent.C results/flight_calibration_research/conditional_mc/reference/extract_parent
