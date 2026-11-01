#!/usr/bin/env bash
set -eo pipefail
cd "$(dirname "$0")/.."
source /nas01/home/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh
native_variant=${1:-light}
case "$native_variant" in
  light) native_source=analysis/tofcalib_supplement_light.C; native_output=results/tofcalib_validation/bin/supplement_light ;;
  full) native_source=analysis/tofcalib_supplement.C; native_output=results/tofcalib_validation/bin/supplement ;;
  probe) native_source=analysis/tofcalib_native_measurement_probe.C; native_output=results/tofcalib_validation/bin/native_measurement_probe ;;
  *) exit 2 ;;
esac
mkdir -p results/tofcalib_validation/bin
"$(root-config --cxx)" -O3 -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0 \
  -D_PGTRACK_ -DUSEGENFIT -D__ROOTSHAREDLIBRARY__ -D_IONL1PRESCALE_ -D_USELAYERRES_ \
  -I../amsd69n -isystem "$AMSWD/include" -isystem "$ROOTSYS/include" -isystem "$AMSWD/GenFit/include" \
  "$native_source" -o "$native_output" \
  -L"$AMSWD/lib/$(root-config --arch)5.34_7" -lntuple_slc6_PG \
  $(root-config --libs) -lMinuit -lTMVA -lNetx -lGeom -lMathMore -lEG -lTreePlayer -lMLP -lXMLIO \
  -L/cvmfs/ams.cern.ch/Offline/CERN/NagLib -lnaggcc64 \
  -L/cvmfs/ams.cern.ch/Offline/CERN/2005.gcc64.44/lib -lmathlib -lgfortran -lm -lgomp
