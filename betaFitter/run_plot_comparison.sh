#!/bin/bash

# Script to run plotBetaComparison.C with parameters from 18 to 26
# Author: GitHub Copilot

echo "Starting to run plotBetaComparison.C with parameters from 18 to 26..."

for i in {18..26}
do
    echo "Running with parameter: $i"
    time root -q "macro/plotBetaComparison.C(\"$i\")"
    echo "Finished parameter: $i"
    echo "----------------------------------------"
done

echo "All runs completed successfully!"