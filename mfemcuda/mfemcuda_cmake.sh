#!/bin/bash

set -x
cd ${STOCKYARD}/mfemcuda
rm -rf build
mkdir build
cd build

if [ "${TACC_SYSTEM}" = "vista" ] ; then
  BLASLIB=${TACC_NVPL_LIB}/libnvpl_blas_lp64_seq.so 
  LAPACKLIB=${TACC_NVPL_LIB}/libnvpl_lapack_lp64_seq.so
else
  # mkl_intel_lp64;mkl_sequential;mkl_core;iomp5;pthread"
  BLASLIB=${TACC_MKL_LIB}/libmkl_sequential.so 
  LAPACKLIB=${TACC_MKL_LIB}/libmkl_intel_lp64.so
fi

CPPSTANDARD=14
export CUDAFLAGS="-allow-unsupported-compiler -Xcompiler -std=gnu++${CPPSTANDARD}"
export NVCC_APPEND_FLAGS="${CUDAFLAGS}"



cmake -Wno-dev \
    -D CMAKE_BUILD_TYPE=Release \
    -D CMAKE_CXX_COMPILER=mpicxx \
    -D CMAKE_CXX_STANDARD=${CPPSTANDARD} \
    \
    -D MFEM_USE_CUDA=YES \
    -D CMAKE_CUDA_FLAGS="${CUDAFLAGS}" \
    \
    -D MFEM_USE_MPI=YES \
    -D MFEM_USE_PETSC=NO \
    -D MFEM_USE_LAPACK=YES \
    -D LAPACK_LIBRARIES=${LAPACKLIB} \
    -D BLAS_LIBRARIES=${BLASLIB} \
    -D MFEM_DEBUG=NO \
    -D HYPRE_DIR=${TACC_HYPRE_DIR} \
    -D MFEM_USE_METIS_5=ON -D METIS_DIR=${TACC_METIS_DIR} \
    -D MFEM_USE_ADIOS2=ON -D ADIOS2_DIR=${TACC_ADIOS2_DIR} \
    ../mfemcuda-4.7

make
