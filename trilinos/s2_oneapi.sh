module unload python2
export TACC_FAMILY_PYTHON_VERSION=3
export TACC_PYTHON_VER=3.7
export TACC_PYTHON_DIR=/opt/apps/intel18/python3/3.7.0
export TACC_PYTHON_INC=/opt/apps/intel18/python3/3.7.0/include
export TACC_PYTHON_BIN=/opt/apps/intel18/python3/3.7.0/bin
export TACC_PYTHON_LIB=/opt/apps/intel18/python3/3.7.0/lib
export TACC_PYTHON_MAN=/opt/apps/intel18/python3/3.7.0/share/man:/opt/apps/intel18/python3/3.7.0/man

module use /scratch/projects/compilers/modulefiles
module load oneapi/22.3.0 oneapi21/impi
export TACC_MKL_DIR=/scratch/projects/compilers/oneapi_2022.3.0/mkl/latest
export TACC_MKL_INC=/scratch/projects/compilers/oneapi_2022.3.0/mkl/latest/include
export TACC_MKL_LIB=/scratch/projects/compilers/oneapi_2022.3.0/mkl/latest/lib/intel46
