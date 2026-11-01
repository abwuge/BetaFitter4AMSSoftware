#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
out=results/flight_calibration_research/mean_loss_likelihood/conditional
mkdir -p "$out"
${CXX:-g++} -O3 -std=c++11 -fopenmp -fPIC -shared analysis/flight_conditional_likelihood.cpp -o "$out/libconditional.so"
sha256sum analysis/flight_conditional_likelihood.cpp "$out/libconditional.so" > "$out/build.sha256"
