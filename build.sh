#!/usr/bin/env bash
set -e

# Setup AMS environment
source /home/ams/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh

export AMSSRC="$AMSWD"

FLAG=""
if [[ "${1:-}" == "debug" ]]; then
    FLAG="debug"
elif [[ "${1:-}" == "clean" ]]; then
    make clean && exit
fi

export USEPRHEION=1
unset USEEVENTORDER 2>/dev/null || true
unset USEHEL1 2>/dev/null || true
unset USENEWL1L9G 2>/dev/null || true
unset USEMCTKRAW 2>/dev/null || true
unset USEONEEV 2>/dev/null || true
unset USEADDTKHIT 2>/dev/null || true
unset USENOLINEARCOR 2>/dev/null || true
unset USECALIB 2>/dev/null || true
export USEHEINNER=1
export USESAVENEG=1
export USEPRL1=1

# Compile AMSVX static library first
echo
echo "build AMSVX..."
echo
current_dir="$(pwd)"
cd "$AMSWD/install/"
if [[ "$FLAG" == "debug" ]]; then
    make debug_static PGTRACK=1 -j$(nproc)
else
    make static PGTRACK=1 -j$(nproc)
fi
cd "$current_dir"

# Build betaFitter
echo
echo "build Pass8..."
echo
export USEPASS7=1

if [[ -n "$FLAG" ]]; then
    make "$FLAG" -j$(nproc)
else
    make all -j$(nproc)
fi
