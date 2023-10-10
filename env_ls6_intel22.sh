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
export VICTOR_WORK=/work2/00434/eijkhout/ls6
export MODULEROOT=${VICTOR_WORK}/modulefiles

export MY_MODULEPATH_ROOT=${MODULEROOT}
module use ${MY_MODULEPATH_ROOT}/Core

intelversion=22.2.0
echo "loading intel ${intelversion}"
module load intel/${intelversion}
module unload intel22/impi
module load impi/22.2.0

module list
