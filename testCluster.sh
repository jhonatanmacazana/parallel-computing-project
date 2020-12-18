#!/bin/bash
#SBATCH -J CPD-test-23:56:22-12/14/20
#SBATCH -p investigacion
#SBATCH -N 1
#SBATCH --tasks-per-node=8

module load gcc/5.5.0

./main.out

module unload gcc/5.5.0
