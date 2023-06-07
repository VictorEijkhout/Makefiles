module reset
module unload intel oneapi gcc impi python2 python3

module use /scratch/projects/compilers/modulefiles

export MODULEROOT=${WORK}/modulefiles
export MY_MODULEPATH_ROOT=${MODULEROOT}

module use ${MY_MODULEPATH_ROOT}/Core

module load intel/23.1.0
module load impi/21.9.0

# module load intel22/python3

module list
