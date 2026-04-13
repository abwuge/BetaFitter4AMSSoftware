#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <base_filename>"
    exit 1
fi

hadd -j32 -f -k results/$1.root results/$1_*.root || { echo "Hadd failed!"; exit 1; }

rm -f results/$1_*.root

echo "Hadd completed!"