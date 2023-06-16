echo "module reset"
module purge
module reset
module list

echo "module unload intel and others"
module unload intel oneapi gcc impi python3
module list

echo "adding experimental and my own module paths"
export MODULEROOT=${WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

intelversion=23.0.0
echo "loading intel ${intelversion}"
module load intel/${intelversion} impi/21.8.0
#module load python3
module list
