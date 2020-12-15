#!/bin/bash

# error flags
set -eu

# change dir to current script
cd "${0%/*}"

# Compilation
make seq

echo "tasks,NRA,time-ms"
echo "sequential"
# Execution
for _i in {2..8}; do
    make run-seq
done

# Cleanup
make clean
