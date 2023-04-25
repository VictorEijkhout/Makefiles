module purge
module reset
module unload intel oneapi impi python2 python3

module use /scratch/projects/compilers/modulefiles

export MODULEROOT=${WORK}/modulefiles
export MY_MODULEPATH_ROOT=${MODULEROOT}

module use ${MY_MODULEPATH_ROOT}/Core

intelversion=22.3.0
module load oneapi/${intelversion}
module load impi/21.4.0

# module load intel22/python3

module list
