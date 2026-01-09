##
## load OneAPI intel compiler
##

echo "==== module reset"
module purge
module reset 2>/dev/null
module load TACC

echo "==== module unload intel and others"
module unload intel oneapi gcc impi mvapich2 python3 python2 2>/dev/null

echo "==== adding experimental and my own module paths"

export VICTOR_WORK=/work/00434/eijkhout/ls6
export VICTOR_MODULEPATH_ROOT=${VICTOR_WORK}/modulefiles
module use ${VICTOR_MODULEPATH_ROOT}/Core

intelversion=25.3
echo "==== loading intel ${intelversion}"
module use /scratch/projects/compilers/modulefiles
module load intel/${intelversion}
module load impi/21.17

##
## Python
##

module load python/3.12

# # get up to date c++ library, which is eff'ed up in the python module
# ## export LD_LIBRARY_PATH=/scratch/tacc/apps/gcc/13.2.0/lib64:${LD_LIBRARY_PATH}

# export TACC_INTEL_LIB=/scratch/projects/compilers/intel24.1/oneapi/2024.1/lib
# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${TACC_INTEL_LIB}

module -t list
