#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <base_filename>"
    exit 1
fi

hadd -j64 -f -k results/$1.root results/$1_*.root

rm -f results/$1_*.root

echo "Hadd completed!"