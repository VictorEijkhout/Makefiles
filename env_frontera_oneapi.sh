echo "module reset"
module purge
module reset
# module list

echo "module unload intel and others"
module unload intel oneapi gcc impi python3
# module list

echo "adding experimental and my own module paths"
module use /scratch1/projects/compilers/modulefiles
export MODULEROOT=${WORK}/modulefiles
export MY_MODULEPATH_ROOT=${MODULEROOT}
module use ${MY_MODULEPATH_ROOT}/Core

intelversion=23.1.0
echo "loading intel ${intelversion}"
module load oneapi/${intelversion}
module load impi/21.9.0
#module load oneapi21/python3
module list
