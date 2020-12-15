#!/bin/bash

# error flags
set -eu

# change dir to current script
cd "${0%/*}"

# Compilation
make par

echo "tasks,N,time-tqli,time-all"
# Execution
for _i in {1..5}; do          # iterations
    make run-par >outputx.csv # code displays taskas, currentN and time
done

# Cleanup
make clean
