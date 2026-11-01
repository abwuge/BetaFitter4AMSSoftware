#!/usr/bin/env bash
set -eo pipefail
cd /nas01/home/hxwu/AMSSoft/betaEstimator
source /nas01/home/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh
"$(root-config --cxx)" -O3 -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0 \
 -D_PGTRACK_ -DUSEGENFIT -D__ROOTSHAREDLIBRARY__ -D_IONL1PRESCALE_ -D_USELAYERRES_ \
 -I../amsd69n -isystem "$AMSWD/include" -isystem "$ROOTSYS/include" -isystem "$AMSWD/GenFit/include" \
 results/terminal_ceiling_20261015/native_readout.C -o results/terminal_ceiling_20261015/extract_native_readout \
 -L"$AMSWD/lib/$(root-config --arch)5.34_7" -lntuple_slc6_PG \
 $(root-config --libs) -lMinuit -lTMVA -lNetx -lGeom -lMathMore -lEG -lTreePlayer -lMLP -lXMLIO \
 -L/cvmfs/ams.cern.ch/Offline/CERN/NagLib -lnaggcc64 \
 -L/cvmfs/ams.cern.ch/Offline/CERN/2005.gcc64.44/lib -lmathlib -lgfortran -lm -lgomp
results/terminal_ceiling_20261015/extract_native_readout \
 /lustre03/data/eos/ams/MC/AMS02/2022/C.B1308/c12.pl1.l1.612000.6_05/82334649.00000001.root \
 results/terminal_ceiling_20261015/native_parts/048.txt \
 results/terminal_ceiling_20261015/native_readout_048.root
