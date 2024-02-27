#!/bin/bash

# Function to execute the main program with specified parameters
execute_main() {
    ./build/main --e 2 -q 4 -b 1000000 -s $1 --graphfile $2 --grid-length $3
}

# Loop for n=20, n=50, and n=70
for n in 20 50 70; do
    # Loop for grid-length 150, 112, 75, and 37
    for grid_length in 150 112 75 37; do
        # Loop for different seeds
        for seed in {100000..100009}; do
            execute_main $seed "data/graph/graph_n${n}_s${seed}.csv" $grid_length
        done
    done
done
