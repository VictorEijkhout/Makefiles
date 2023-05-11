module reset
module unload intel gcc impi python2 python3

module use /scratch/projects/compilers/modulefiles

export MODULEROOT=${WORK}/modulefiles
export MY_MODULEPATH_ROOT=${MODULEROOT}

module use ${MY_MODULEPATH_ROOT}/Core

module load intel/19.1.1

module load impi/19.0.9

