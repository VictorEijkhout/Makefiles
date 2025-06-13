##
## load OneAPI intel compiler
##

echo "module reset"
module purge
module reset
# module list
module load TACC

echo "module unload intel and others"
module unload intel oneapi gcc nvidia impi mvapich2 openmpi python3 python2 2>/dev/null

export VICTOR_WORK=/work2/00434/eijkhout/stampede3
export MODULEROOT=${VICTOR_WORK}/modulefiles

export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

nvversion=25.3
echo "loading nvidia ${nvversion}"
module load nvidia/${nvversion}

module load openmpi

module -t list 2>&1 | sort | awk '{v=v" "$0} END {print v}'

