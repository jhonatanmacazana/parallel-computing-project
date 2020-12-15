#!/bin/bash

# error flags
set -eu

# change dir to current script
cd "${0%/*}"

# send Cluster
tar -zcf temp.tar.gz Makefile testLocal.sh src/parallel/{main.cpp,utils.hpp}
scp temp.tar.gz cpd:~/proj

rm temp.tar.gz

# # Compilation
# make par

# echo "tasks,N,time-tqli,time-all"
# # Execution
# for _i in {1..100}; do # iterations
#     make run-par       # code displays taskas, currentN and time
# done

# # Cleanup
# make clean
