module purge
module reset
module unload intel oneapi gcc impi python2 python3 2>/dev/null

## module use /scratch/projects/compilers/modulefiles

export VICTOR_WORK=/work2/00434/eijkhout/stampede2
export MODULEROOT=${VICTOR_WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}

module use ${VICTOR_MODULEPATH_ROOT}/Core

gccversion=12.2.0
module load gcc/${gccversion}
export TACC_CC=gcc
export TACC_CXX=g++
export TACC_FC=gfortran
module load impi/19.0.9

module list
