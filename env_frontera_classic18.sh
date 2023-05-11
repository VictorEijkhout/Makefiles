echo "module reset"
module purge
module reset
module list

echo "module unload intel and others"
module unload intel oneapi gcc impi python3
module list

echo "adding experimental and my own module paths"
## THIS CONTAINS CRAP module use /scratch1/projects/compilers/modulefiles
export MODULEROOT=${WORK}/modulefiles
export MY_MODULEPATH_ROOT=${MODULEROOT}
module use ${MY_MODULEPATH_ROOT}/Core

intelversion=18.0.5
echo "loading intel ${intelversion}"
module load intel/${intelversion} impi/${intelversion}
module load python3
module list
