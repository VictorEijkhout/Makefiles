echo "module reset"
module reset
module list

echo "module unload intel and others"
module unload intel oneapi gcc impi python3
module list

echo "adding experimental and my own module paths"
module use /scratch1/projects/compilers/modulefiles
export MODULEROOT=${WORK}/modulefiles
export MY_MODULEPATH_ROOT=${MODULEROOT}
module use ${MY_MODULEPATH_ROOT}/Core

echo "loading intel 22.3"
module load intel/22.3.0 impi/22.3.0
#module load oneapi21/python3
module list
