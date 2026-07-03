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

intelversion=24.0
echo "Available intel/${intelversion}:"
module -t avail intel/${intelversion}
echo " .. loading intel ${intelversion}"

module load intel/${intelversion}
module -t show intel

export TACC_CC=icx
export TACC_CXX=icpx
export TACC_FC=ifx
module load impi/21.11

echo "Module path:"
echo $MODULEPATH | tr ':' '\n'
echo "Loaded:"
module -t list
echo "Intel:"
module -t show intel
