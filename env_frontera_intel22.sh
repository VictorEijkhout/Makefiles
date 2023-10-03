echo "module reset"
module reset

echo "module unload intel and others"
module unload intel oneapi gcc impi python3 2>/dev/null

echo "adding experimental and my own module paths"
module use /scratch1/projects/compilers/modulefiles
export VICTOR_WORK=/work2/00434/eijkhout/frontera
export MODULEROOT=${VICTOR_WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

echo "loading intel 22.3"
module load intel/22.3.0 impi/22.3.0

module load oneapi21/python3

module list
