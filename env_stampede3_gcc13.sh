##
## load OneAPI intel compiler
##

echo "module reset"
module purge
module reset
# module list
module load TACC

echo "module unload intel and others"
module unload intel oneapi gcc impi mvapich2 python3 python2 2>/dev/null
# module list

export VICTOR_WORK=/work2/00434/eijkhout/stampede3
export MODULEROOT=${VICTOR_WORK}/modulefiles

export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

gccversion=13.2.0
impiversion=21.9
echo "loading gcc ${gccversion}"
module load gcc/${gccversion}
module load mkl
module unload impi
module load impi/${impiversion}
# /21.9 or 21.11?
export PATH=${HOME}/bin/gcc15:${PATH}

module -t list | sort

export PATH=${STOCKYARD}/MrPackMod:${PATH}
export PYTHONPATH=${STOCKYARD}:${PYTHONPATH}
