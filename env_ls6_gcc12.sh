##
## load OneAPI intel compiler
##

echo "module reset"
module purge
module reset
# module list
module load TACC

echo "module unload intel and others"
module unload intel oneapi gcc impi mvapich2 python3 python2 2>/dev/null 2>/dev/null
# module list

echo "adding experimental and my own module paths"
## module use /scratch/projects/compilers/modulefiles
export VICTOR_WORK=/work2/00434/eijkhout/ls6
export MODULEROOT=${VICTOR_WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

gccversion=12.2.0
echo "loading gcc ${gccversion}"
module load gcc/${gccversion}
module load impi/19.0.9

module list
