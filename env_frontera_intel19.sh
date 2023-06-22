echo "module reset"
module reset
module list

echo "module unload intel and others"
module unload intel oneapi gcc impi python3 2>/dev/null
module list

echo "adding experimental and my own module paths"
export MODULEROOT=${WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

intelversion=19.1.1
echo "loading intel ${intelversion}"
module load intel/${intelversion} impi/19.0.9
module load python3

module list
