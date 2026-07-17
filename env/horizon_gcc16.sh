##
## load GCC 16
##

echo "module reset"
module -t purge
module -t reset
# module -t list
module -t load TACC

echo "module unload intel and others"
module -t unload intel oneapi gcc impi mvapich2 python3 python2 2>/dev/null

export VICTOR_WORK=${WORK}
export MODULEROOT=${VICTOR_WORK}/modulefiles

export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

gccversion=16.1.0
echo "loading gcc ${gccversion}"
module -t load gcc/${gccversion}
module -t load nvpl
module -t load openmpi
# ??? export PATH=${HOME}/bin/gcc15:${PATH}

# append for libisl
# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/apps/gcc/15.1.0/lib

module -t list 2>&1 | sort | awk '{v=v" "$0} END {print v}'

export PATH=${STOCKYARD}/MrPackMod:${PATH}
export PYTHONPATH=${STOCKYARD}:${PYTHONPATH}
