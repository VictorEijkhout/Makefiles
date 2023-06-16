module reset
module unload intel gcc impi python3 python2

module use /scratch/projects/compilers/modulefiles

export MODULEROOT=${WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}

module use ${VICTOR_MODULEPATH_ROOT}/Core

module load intel/22.2.0

# TACC provides: intel22/impi/22.2.0
module load impi/22.2.0

module load intel22/python3

module list
