#!/bin/bash

# error flags
set -eu

# change dir to current script
cd "${0%/*}"

# compress files
tar -zcf temp.tar.gz Makefile testLocal.sh src/{parallel,sequential}/{main.cpp,utils.hpp}

# send to  Cluster
scp temp.tar.gz cpd:~/proj

# cleanup
rm temp.tar.gz
