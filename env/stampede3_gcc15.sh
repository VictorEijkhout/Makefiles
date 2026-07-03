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

export VICTOR_WORK=/work2/00434/eijkhout/stampede3
export MODULEROOT=${VICTOR_WORK}/modulefiles

export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

gccversion=15.1.0
echo "loading gcc ${gccversion}"
module load gcc/${gccversion}
module load mkl
module load impi
export PATH=${HOME}/bin/gcc15:${PATH}

# append for libisl
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/apps/gcc/15.1.0/lib

module -t list 2>&1 | sort | awk '{v=v" "$0} END {print v}'

export PATH=${STOCKYARD}/MrPackMod:${PATH}
export PYTHONPATH=${STOCKYARD}:${PYTHONPATH}
