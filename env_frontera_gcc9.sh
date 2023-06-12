echo "module reset"
module purge
module reset
module list

echo "module unload intel and others"
module unload intel oneapi impi python3
module list

echo "adding experimental and my own module paths"
export MODULEROOT=${WORK}/modulefiles
export MY_MODULEPATH_ROOT=${MODULEROOT}
module use ${MY_MODULEPATH_ROOT}/Core

gnuversion=9.1.0
echo "loading gnu ${gnuversion}"
module load gcc/${gnuversion} impi/21.9.0
#module load python3
module list
