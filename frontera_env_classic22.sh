module reset
module unload intel impi python3

module use /scratch1/projects/compilers/modulefiles

export MODULEROOT=${WORK}/modulefiles
export MY_MODULEPATH_ROOT=${MODULEROOT}

module use ${MY_MODULEPATH_ROOT}/Core

module load intel/22.3.0 impi/22.3.0