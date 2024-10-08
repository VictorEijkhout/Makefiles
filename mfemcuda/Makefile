################################################################
####
#### Makefile for Mfem installation
####
#### only needed for the CHyPS CSA project
####
################################################################

PACKAGE = mfemcuda
URL = https://github.com/mfem/mfem
PACKAGEVERSION = 4.7
## 4.5.2 has breaking changes for CHyPS.
HYPREVERSION = 2.31.0

MODE = mpi
BASICMODULES = mkl adios2 hypre/${HYPREVERSION}-i64 metis 
MODULES = ${BASIC_MODULES} cuda

include ${MAKEINCLUDES}/Make.info

##
## cmake 
##
ifeq "${TACC_SYSTEM}" "vista"
  BLASLIB=${TACC_NVPL_LIB}/libnvpl_blas_lp64_seq.so 
  LAPACKLIB=${TACC_NVPL_LIB}/libnvpl_lapack_lp64_seq.so
else
  # mkl_intel_lp64;mkl_sequential;mkl_core;iomp5;pthread"
  BLASLIB=${TACC_MKL_LIB}/libmkl_sequential.so 
  LAPACKLIB=${TACC_MKL_LIB}/libmkl_intel_lp64.so
endif

USEGPU = YES
CUDAFLAGS = -allow-unsupported-compiler -Xcompiler -std=gnu++14
CMAKEFLAGS = -Wno-dev \
    -D CMAKE_BUILD_TYPE=Release \
    -D MFEM_USE_MPI=YES \
    -D MFEM_USE_CUDA=${USEGPU} \
    -D CMAKE_CUDA_FLAGS="-std=c++14" \
    -D MFEM_USE_PETSC=NO \
    -D MFEM_USE_LAPACK=YES \
    -D LAPACK_LIBRARIES=${LAPACKLIB} \
    -D BLAS_LIBRARIES=${BLASLIB} \
    -D MFEM_DEBUG=NO \
    -D HYPRE_DIR=${TACC_HYPRE_DIR} \
    -D MFEM_USE_METIS_5=ON -D METIS_DIR=${TACC_METIS_DIR} \
    -D MFEM_USE_ADIOS2=ON -D ADIOS2_DIR=${TACC_ADIOS2_DIR}
CMAKEPREFIXPATHSET=1

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : gpu version"
.PHONY: default_install gpu
default_install gpu :
	@export CUDAFLAGS="${CUDAFLAGS}" \
	 && export NVCC_APPEND_FLAGS="${CUDAFLAGS}" \
	 && make --no-print-directory configure build 

TGZURL = https://github.com/mfem/mfem/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = git@github.com:mfem/mfem.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean
