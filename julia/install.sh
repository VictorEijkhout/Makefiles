#!/bin/bash

#SBATCH -J myjob           # Job name
#SBATCH -o myjob.o%j       # Name of stdout output file
#SBATCH -e myjob.e%j       # Name of stderr error file
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 10:0:0
#SBATCH -A A-ccsc       # Allocation name (req'd if you have more than 1)

PACKAGEROOT=$WORK  make default_install PACKAGEVERSION=git JCOUNT=4
