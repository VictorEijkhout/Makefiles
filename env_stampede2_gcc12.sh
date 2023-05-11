module purge
module reset
module unload intel oneapi gcc impi python2 python3

## module use /scratch/projects/compilers/modulefiles

export MODULEROOT=${WORK}/modulefiles
export MY_MODULEPATH_ROOT=${MODULEROOT}

module use ${MY_MODULEPATH_ROOT}/Core

gccversion=12.2.0
module load gcc/${gccversion}
module load impi/19.0.5

module list
