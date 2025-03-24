#!/bin/bash

# Default maximum number of parallel processes
Z=${1:-8}
fitOption=${2:-0}
energyLossScale=${3:-1.0}
MAX_PROCS=${4:-100}

# Array to store child PIDs
declare -a CHILD_PIDS

# Function to cleanup child processes
cleanup() {
    echo -e "\nReceived termination signal. Cleaning up..."
    # Kill all child processes
    for pid in "${CHILD_PIDS[@]}"; do
        if kill -0 $pid 2>/dev/null; then
            kill $pid 2>/dev/null
        fi
    done
    exit 1
}

# Set up signal handlers
trap cleanup SIGINT SIGTERM

# Get the script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
INPUT_LIST="${SCRIPT_DIR}/input_Z${Z}.list"
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
    ("${SCRIPT_DIR}/run.csh" "$input_file" "$output_file" "$fitOption" "$energyLossScale" > "$log_file" 2>&1) &
    
    # Store the PID of the background process
    CHILD_PIDS+=($!)
    
    echo "Started processing: $input_file"
done < "$INPUT_LIST"

# Wait for all background processes to complete
wait

read -t 0 -p "Hadd all root files? [Y/n] " response
response=${response:-Y}
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]; then
    ./hadd.sh $RUN_NUM
else
    echo "Skipping hadd operation."
fi

# Add information to README.md
README_FILE="${RESULTS_DIR}/README.md"
TIMESTAMP=$(date "+%Y-%m-%d %H:%M:%S")
HADD_FILE="${RUN_NUM}.root"

# Create README.md if it doesn't exist
if [ ! -f "$README_FILE" ]; then
    echo "# Running Parameters Log" > "$README_FILE"
fi

# Append run information to README.md
echo "[${TIMESTAMP}] FILE = ${HADD_FILE}: Z = ${Z}, fitOption = ${fitOption}, energyLossScale = ${energyLossScale}, MAX_PROCS = ${MAX_PROCS}  " >> "$README_FILE"
echo "Run information has been added to ${README_FILE}."

echo "All jobs completed!"