echo "module reset"
module reset
## module list

echo "module unload intel and others"
module unload intel oneapi gcc impi python2 python3 2>/dev/null

echo "adding experimental and my own module paths"
export VICTOR_WORK=/work2/00434/eijkhout/frontera
export MODULEROOT=${VICTOR_WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

gnuversion=9.1.0
echo "loading gnu ${gnuversion}"
module load gcc/${gnuversion} 
export TACC_CC=gcc
export TACC_CXX=g++
export TACC_FC=gfortran
module load impi/19.0.9

#module load python3
module list
