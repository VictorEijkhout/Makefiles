##
## load OneAPI intel compiler
##

echo "module reset"
module purge
module reset
# module list

echo "module unload intel and others"
module unload intel oneapi gcc impi mvapich2 python3 python2 2>/dev/null
module use  /scratch1/projects/compilers/pvcmodules
module use  /scratch1/projects/compilers/modulefiles

echo "adding experimental and my own module paths"
## ONEAPI IS OFFICIAL as Intel module
export VICTOR_WORK=/work2/00434/eijkhout/frontera
export MODULEROOT=${VICTOR_WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

echo "loading intel ${intelversion}"
module load oneapi/24
module load oneapi22/impi/22.3.0

module load oneapi21/onepython3/3.7.0

# export PATH=/opt/apps/intel19/python3/${pvv}/bin:${PATH}
# export LD_LIBRARY_PATH=/opt/apps/intel19/python3/${pvv}/lib:${LD_LIBRARY_PATH}
# export PYTHONPATH=/opt/apps/intel19/impi19_0/python3/${pvv}/lib/python${pv}/site-packages:${PYTHONPATH}

module list
