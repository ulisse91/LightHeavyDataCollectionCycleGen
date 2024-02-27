#!/bin/bash

set -e

timestamp=$(date "+%Y%m%d%H%M")
main="./build/main"
pathoutput="data/output/"
pathbkp="data/bkp/"
baseseed=100000

function check_command {
    if [ $2 -eq 0 ]; then
        printf "%s [OK]\n" "$1"
    else
        printf "%s [FAILED]\n" "$1"
        cleanup "$2"
    fi
}

function check_error {
    if [ $2 -ne 0 ]; then
        printf "%s [ERROR]\n" "$1"
        cleanup "$2"
    fi
}

cleanup() {
    local exit_code=$1
    printf "Error occurred during execution (Exit Code: %d)\n" "$exit_code"
    printf "Cleaning up...\n"

    # Kill all the background processes
    for pid in ${pids[@]}; do
        kill $pid
    done

    rm -rf "data/output"
    exit "$exit_code"
}

trap 'cleanup $?' ERR INT
trap 'cleanup $?' SIGINT

if [[ "$1" -eq 1 ]]; then
    printf "Starting compilation ...\n\n"
    make main
    printf "\n"
    check_command "Compilation" $?
fi

mkdir -p "$pathbkp"
mkdir -p "$pathoutput"

printf "\nStarting simulations ... \n"
i=0
drone=4
for seed in {0..9}; do
    for nodes in {20,50,70}; do
        for budget in {750000,1000000,1500000}; do
            "$main" --e 1 -b "$budget" -q "$drone" -n "$nodes" -s $((baseseed + seed)) >> "${pathoutput}results-b${budget}-d${drone}-n${nodes}.txt" &
            pids[$budget]=$!  # Store the process ID of each background process
        done

        # Wait for all the background processes to complete
        for pid in ${pids[@]}; do
            wait $pid
            check_error "Execution" $?

            ((i += 1))
            progress=$((i * 100 / (3 * 3 * 10)))  # Calculate the percentage progress
            printf "Progress: %d%%\r" "$progress"
        done

        unset pids  # Clear the array of process IDs

        # Sequential code for running the experiments:
        # for seed in {0..9}; do
        #     "$main" --e 1 -b "$budget" -q "$drone" -n "$nodes" -s $((baseseed + seed)) >> "${pathoutput}results-b${budget}-d${drone}-n${nodes}.txt"
        #     check_error "Execution" $?
        #     ((i += 1))
        #     progress=$((i * 100 / (3 * 3 * 10)))  # Calculate the percentage progress
        #     printf "Progress: %d%%\r" "$progress"
        # done
    done
done


printf "\nSimulation terminated\n"

printf "Starting scripts ...\n"
python3 scripts/filter_files_cycles.py data/output/ data/output/ A B C E F G H I L M N O P Q R S T U V 0


python3 scripts/read_output.py

printf "Starting compression files ... "
tar -czf "${pathbkp}${timestamp}.tar.gz" "${pathoutput}"

printf "DONE\n"