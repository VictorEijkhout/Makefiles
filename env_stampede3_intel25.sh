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

intelversion=25.0
echo "loading intel ${intelversion}"

export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
## 25 is on scratch
module use /scratch/projects/compilers/modulefiles
module use ${VICTOR_MODULEPATH_ROOT}/Core

## echo $MODULEPATH | tr ':' '\n'

module load intel/${intelversion}
# export TACC_CC=icx
# export TACC_CXX=icpx
# export TACC_FC=ifx
module load impi

echo "Loaded:"
module -t list
