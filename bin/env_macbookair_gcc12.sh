##
## load OneAPI intel compiler
##

echo "module reset"
module purge
module reset
# module list
# module load TACC

echo "module unload intel and others"
module unload intel oneapi gcc impi mvapich2 python3 python2 2>/dev/null 2>/dev/null
# module list

echo "adding experimental and my own module paths"
export MODULEROOT=${HOME}/modulefiles
export MY_MODULEPATH_ROOT=${HOME}/modulefiles
module use ${MY_MODULEPATH_ROOT}/Core

gccversion=12
echo "loading gcc ${gccversion}"
module load gcc/${gccversion}
module load mpich

module list
