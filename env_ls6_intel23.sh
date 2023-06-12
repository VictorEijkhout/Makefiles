##
## load OneAPI intel compiler
##

echo "module reset"
module purge
module reset
# module list
module load TACC

echo "module unload intel and others"
module unload intel oneapi gcc impi mvapich2 python3 python2
# module list

echo "adding experimental and my own module paths"
## ONEAPI IS not yet OFFICIAL as INTEL
module use /scratch/projects/compilers/modulefiles
export MODULEROOT=${WORK}/modulefiles
export MY_MODULEPATH_ROOT=${MODULEROOT}
module use ${MY_MODULEPATH_ROOT}/Core

intelversion=23.1.0
echo "loading intel ${intelversion}"
module load intel/${intelversion}
module load impi/21.9.0

module list