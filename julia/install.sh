#!/bin/bash

#SBATCH -J myjob           # Job name
#SBATCH -o myjob.o%j       # Name of stdout output file
#SBATCH -e myjob.e%j       # Name of stderr error file
#SBATCH -p skx
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 1              # Total # of mpi tasks
#SBATCH -t 10:0:0
#SBATCH -A A-ccsc       # Allocation name (req'd if you have more than 1)

PACKAGEROOT=$WORK  make default_install PACKAGEVERSION=git JCOUNT=4
