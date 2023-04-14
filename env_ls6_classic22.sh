module reset
module unload intel impi python3

module use /scratch/projects/compilers/modulefiles

export MODULEROOT=${WORK}/modulefiles
export MY_MODULEPATH_ROOT=${MODULEROOT}

module use ${MY_MODULEPATH_ROOT}/Core

module load intel/22.2.0

# TACC provides: intel22/impi/22.2.0
module load impi/22.2.0

module load intel22/python3
