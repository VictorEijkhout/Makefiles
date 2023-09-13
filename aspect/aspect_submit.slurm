#!/bin/bash
#----------------------------------------------------
#SBATCH -J myjob           # Job name
#SBATCH -o myjob.o%j       # Name of stdout output file
#SBATCH -e myjob.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 06:30:00        # Run time (hh:mm:ss)
#SBATCH -A A-ccsc
# #SBATCH --mail-type=all
#SBATCH --mail-user=eijkhout@tacc.utexas.edu

# Launch serial code...
./aspect_ladder -j 40 11

