#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
out=results/flight_calibration_research/continuation_20261013/periodic_adc/conditional_noise
mkdir -p "$out"
if [[ -e "$out/libconditionalnoise.so" ]]; then
  echo 'Refusing to replace an existing conditional-noise binary; inspect its freeze first.' >&2
  exit 1
fi
g++ -O3 -std=c++11 -fPIC -shared -fopenmp analysis/conditional_time_noise.cpp -o "$out/libconditionalnoise.so"
sha256sum "$out/libconditionalnoise.so"
