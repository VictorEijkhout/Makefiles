module reset
module unload intel gcc impi python3 python2 2>/dev/null

module use /scratch/projects/compilers/modulefiles

export VICTOR_WORK=/work/00434/eijkhout/ls6
export MODULEROOT=${VICTOR_WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}

module use ${VICTOR_MODULEPATH_ROOT}/Core

module load intel/22.2.0

# TACC provides: intel22/impi/22.2.0
module load impi/22.2.0

module load intel22/python3

module list
