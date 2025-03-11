#!/bin/bash

# Default maximum number of parallel processes
MAX_PROCS=${1:-32}

# Get the script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
INPUT_LIST="${SCRIPT_DIR}/input.list"
RESULTS_DIR="${SCRIPT_DIR}/results"
LOGS_DIR="${SCRIPT_DIR}/logs"

# Create directories if they don't exist
mkdir -p "$RESULTS_DIR"
mkdir -p "$LOGS_DIR"

# Function to get number of running jobs
get_running_jobs() {
    jobs -p | wc -l
}

# Process each input file
while IFS= read -r input_file; do
    # Skip empty lines and comments
    [[ -z "$input_file" || "${input_file:0:1}" == "#" ]] && continue
    
    # Wait if we have reached the maximum number of processes
    while [ $(get_running_jobs) -ge $MAX_PROCS ]; do
        sleep 1
    done
    
    # Generate output filename based on input filename
    filename=$(basename "$input_file")
    output_file="${RESULTS_DIR}/${filename}"
    log_file="${LOGS_DIR}/${filename%.root}.log"
    
    # Run the process in background
    ("${SCRIPT_DIR}/run.csh" "$input_file" "$output_file" > "$log_file" 2>&1) &
    
    echo "Started processing: $input_file"
done < "$INPUT_LIST"

# Wait for all background processes to complete
wait

echo "All jobs completed!"