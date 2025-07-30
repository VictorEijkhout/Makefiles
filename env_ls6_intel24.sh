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

export VICTOR_WORK=/work2/00434/eijkhout/ls6
export VICTOR_MODULEPATH_ROOT=${VICTOR_WORK}/modulefiles
module use ${VICTOR_MODULEPATH_ROOT}/Core

intelversion=24.1
echo "==== loading intel ${intelversion}"
module load intel/${intelversion}
module load impi/21.12

##
## Python
##

# get our python somewhere
export TACC_PYTHON_DIR=/opt/apps/gcc11_2/python3/3.9.7
# export TACC_PYTHON_DIR=/opt/apps/intel19/python3/3.9.7

export TACC_PYTHON3_DIR=${TACC_PYTHON_DIR}
export PATH=${TACC_PYTHON_DIR}/bin:${PATH}
export LD_LIBRARY_PATH=${TACC_PYTHON_DIR}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${TACC_PYTHON_DIR}/lib:${LD_LIBRARY_PATH}

# ## export PYTHONPATH=${PYTHONPATH}:/opt/apps/intel19/impi19_0/python3/3.9.7/lib/python3.9/site-packages

# # get up to date c++ library, which is eff'ed up in the python module
# ## export LD_LIBRARY_PATH=/scratch/tacc/apps/gcc/13.2.0/lib64:${LD_LIBRARY_PATH}

# export TACC_INTEL_LIB=/scratch/projects/compilers/intel24.1/oneapi/2024.1/lib
# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${TACC_INTEL_LIB}

module -t list
