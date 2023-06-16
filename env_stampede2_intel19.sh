module reset
module unload intel oneapi gcc impi python2 python3 2>/dev/null

## module use /scratch/projects/compilers/modulefiles

export MODULEROOT=${WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}

module use ${VICTOR_MODULEPATH_ROOT}/Core

module load intel/19.1.1
module load impi/19.0.9

# module load intel22/python3

module list
