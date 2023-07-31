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

echo "adding experimental and my own module paths"
## ONEAPI IS not yet OFFICIAL as INTEL
module use /scratch/projects/compilers/modulefiles
export MODULEROOT=${WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

intelversion=19.1.1
echo "loading intel ${intelversion}"
module load intel/${intelversion}
module load impi/19.0.9
module load python3

module list