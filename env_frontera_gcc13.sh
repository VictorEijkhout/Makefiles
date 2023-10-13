echo "module reset"
module purge
module reset
## module list

echo "module unload intel and others"
module unload intel oneapi impi python3 2>/dev/null
## module list

echo "adding experimental and my own module paths"
export MODULEROOT=/work2/00434/eijkhout/frontera/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

gnuversion=13.2.0
echo "loading gnu ${gnuversion}"
module load gcc/${gnuversion} 
export TACC_CC=gcc
export TACC_CXX=g++
export TACC_FC=gfortran

module load impi/21.9.0

module load mkl

#module load python3
module list
