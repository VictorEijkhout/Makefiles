echo "module reset"
export MODULEPATH=$( echo $MODULEPATH | sed -e 's/:/\n/' | grep -v eijkhout | assemblepath )
module reset
module purge 
module load TACC

echo "module unload intel and others"
module unload intel oneapi gcc impi python2 python3 2>/dev/null

echo "adding experimental and my own module paths"
export VICTOR_WORK=/work2/00434/eijkhout/frontera
export MODULEROOT=${VICTOR_WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

intelversion=19.1.1
echo "loading intel ${intelversion}"
module load intel/${intelversion} impi/19.0.9
export TACC_CC=icc
export TACC_CXX=icpc
export TACC_FC=ifort
module load python3

module list
