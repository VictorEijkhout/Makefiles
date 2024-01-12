##
## load OneAPI intel compiler
##

echo "module reset"
module purge
module reset
# module list
module load TACC

echo "module unload intel and others"
module unload intel oneapi gcc impi mvapich2 python3 python2 2>/dev/null
# module list

export VICTOR_WORK=/work2/00434/eijkhout/stampede3
export MODULEROOT=${VICTOR_WORK}/modulefiles

export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

intelversion=24.0
echo "loading intel ${intelversion}"
module load intel/${intelversion}
export TACC_CC=icx
export TACC_CXX=icpx
export TACC_FC=ifx
module load impi/21.11

module list
