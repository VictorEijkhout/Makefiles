#!/bin/bash

#SBATCH -J installjob           # Job name
#SBATCH -o installjob.o%j       # Name of stdout output file
#SBATCH -e installjob.o%j       # Name of stderr error file
#SBATCH -p skx-dev             # Queue (partition) name
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 02:0:0
##SBATCH --mail-user=myname@myschool.edu
##SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A A-ccsc       # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...

cd cxxopts
make default_install JCOUNT=20
