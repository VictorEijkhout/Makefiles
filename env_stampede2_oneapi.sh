module purge
module reset
module unload intel oneapi gcc impi python2 python3 2>/dev/null

module use /scratch/projects/compilers/modulefiles

export VICTOR_WORK=/work2/00434/eijkhout/stampede2
export MODULEROOT=${VICTOR_WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}

module use ${VICTOR_MODULEPATH_ROOT}/Core

intelversion=22.3.0
module load oneapi/${intelversion}
export TACC_CC=icx
export TACC_CC=icpx
export TACC_FC=ifx
module load impi/21.4.0

# module load intel22/python3

module list
