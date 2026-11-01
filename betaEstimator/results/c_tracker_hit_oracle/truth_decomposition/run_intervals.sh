#!/usr/bin/env bash
set -eo pipefail
source /nas01/home/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh >/dev/null
cd /nas01/home/hxwu/AMSSoft/betaEstimator/results/c_tracker_hit_oracle/truth_decomposition
g++ -O2 -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0 -D_PGTRACK_ -DUSEGENFIT -D__ROOTSHAREDLIBRARY__ -I"$AMSWD/include" -I"$ROOTSYS/include" -I"$AMSWD/GenFit/include" transport_intervals.C -o transport_intervals -L"$AMSWD/lib/linuxx8664gcc5.34_7" -lntuple_slc6_PG $(root-config --libs) -lMinuit -lTMVA -lNetx -lGeom -lMathMore -lEG -lTreePlayer -lMLP -lXMLIO -L/cvmfs/ams.cern.ch/Offline/CERN/NagLib -lnaggcc64 -L/cvmfs/ams.cern.ch/Offline/CERN/2005.gcc64.44/lib -lmathlib -lgfortran -lm -lgomp
./transport_intervals
