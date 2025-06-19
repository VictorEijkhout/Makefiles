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
export VICTOR_MODULEPATH_ROOT=${VICTOR_WORK}/modulefiles
module use ${VICTOR_MODULEPATH_ROOT}/Core

intelversion=23.1
echo "loading intel ${intelversion}"

module load intel/${intelversion}
export TACC_CC=icx
export TACC_CXX=icpx
export TACC_FC=ifx
module load impi/21.13

echo "Loaded:"
module -t list
