##
## load OneAPI intel compiler
##

echo "==== module reset"
module purge
module reset 2>/dev/null
# module list
module load TACC

echo "==== module unload intel and others"
module unload intel oneapi gcc impi mvapich2 python3 python2 2>/dev/null
# module list

echo "==== adding experimental and my own module paths"

export VICTOR_WORK=/work2/00434/eijkhout/ls6
export MODULEROOT=${VICTOR_WORK}/modulefiles

export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

intelversion=24.1
echo "==== loading intel ${intelversion}"
module load intel/${intelversion}
# export TACC_CC=icx
# export TACC_CXX=icpx
# export TACC_FC=ifx
module load impi/21.12

## let's see if this python is compatible enough
## module load oneapi22/python3/3.9.10
module use /scratch/projects/compilers/modulefiles
module load intel22/python3

module -t list
