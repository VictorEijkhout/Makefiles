##
## load OneAPI intel compiler
##

echo "module reset"
module purge
module reset
# module list
module load TACC

echo "module unload intel and others"
module unload intel oneapi gcc impi mvapich2 python3 python2 2>/dev/null

echo "adding experimental and my own module paths"
## ONEAPI IS OFFICIAL as INTEL module use /scratch1/projects/compilers/modulefiles
export VICTOR_WORK=/work2/00434/eijkhout/frontera
export MODULEROOT=${VICTOR_WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

intelversion=23.1.0
echo "loading intel ${intelversion}"
module load intel/${intelversion}
module load impi/21.9.0

module list

echo "add python3 from intel21"
export PATH=/scratch1/projects/compilers/oneapi_2021.4.0.3422/intelpython/python3.7/bin/:$PATH
export PYTHONPATH=/opt/apps/intel19/impi19_0/python3/3.7.0/lib/python3.7/site-packages:$PYTHONPATH
