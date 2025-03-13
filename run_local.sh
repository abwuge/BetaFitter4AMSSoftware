#!/bin/bash

# Default maximum number of parallel processes
MAX_PROCS=${1:-64}

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

get_next_run_num() {
    local i=0
    while [[ -e "${RESULTS_DIR}/${i}.root" ]]; do
        i=$((i + 1))
    done
    echo "$i"
}
RUN_NUM=$(get_next_run_num)

# Process each input file
counter=0
while IFS= read -r input_file; do
    # Skip empty lines and comments
    [[ -z "$input_file" || "${input_file:0:1}" == "#" ]] && continue
    
    # Wait if we have reached the maximum number of processes
    while [ $(get_running_jobs) -ge $MAX_PROCS ]; do
        sleep 1
    done
    
    # Generate output filename
    basename=${RUN_NUM}_${counter}
    output_file="${RESULTS_DIR}/${basename}.root"
    log_file="${LOGS_DIR}/${basename}.log"
    
    ((counter++))
    
    # Run the process in background
    ("${SCRIPT_DIR}/run.csh" "$input_file" "$output_file" > "$log_file" 2>&1) &
    
    echo "Started processing: $input_file"
done < "$INPUT_LIST"

# Wait for all background processes to complete
wait

read -p "hadd all root files? [Y/n] " response
response=${response:-Y}
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]; then
    ./hadd.sh $RUN_NUM
else
    echo "Skipping hadd operation."
fi

echo "All jobs completed!"