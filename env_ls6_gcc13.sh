##
## load OneAPI intel compiler
##

echo "module reset"
module purge
module reset
# module list
module load TACC

echo "module unload intel and others"
module unload intel oneapi gcc impi mvapich2 python3 python2 2>/dev/null 2>/dev/null
# module list

echo "adding experimental and my own module paths"
## module use /scratch/projects/compilers/modulefiles
export VICTOR_WORK=/work2/00434/eijkhout/ls6
export MODULEROOT=${VICTOR_WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

compiler=gcc
compilerversion=13.2.0
echo "loading ${compiler} ${compilerversion}"
module load ${compiler}/${compilerversion}
export TACC_CC=gcc
export TACC_CXX=g++
export TACC_FC=gfortran
module load impi
module load mkl

# module load python3/3.9.7
export TACC_PYTHON_DIR=/opt/apps/gcc11_2/python3/3.9.7
export TACC_PYTHON3_DIR=/opt/apps/gcc11_2/python3/3.9.7
export PATH=${TACC_PYTHON_DIR}/bin:${PATH}
export TACC_INTEL_LIB=/scratch/projects/compilers/intel24.1/oneapi/2024.1/lib
export LD_LIBRARY_PATH=${TACC_PYTHON_DIR}/lib:${LD_LIBRARY_PATH}:${TACC_INTEL_LIB}

module list
