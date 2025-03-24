#!/bin/tcsh

source /nas01/home/hxwu/AMSSoft/amsroot534.csh

if ($# < 2) then
    echo "Error: Need at least 2 parameters"
    echo "Usage: $0 inputFile outputFile [fitOption=-2] [energyLossScale=1.0]"
    exit 1
endif

set inputFile = "$1"
set outputFile = "$2"

# Handle fitOption (3rd parameter with default -2)
if ($# >= 3) then
    set fitOption = "$3"
else
    set fitOption = "0"
endif

# Handle energyLossScale (4th parameter with default 2.0)
if ($# >= 4) then
    set energyLossScale = "$4"
else
    set energyLossScale = "1.0"
endif

./build/bin/run/betaFitter $inputFile $outputFile $fitOption $energyLossScale