################################################################
####
#### Makefile for Mfem installation
####
################################################################

PACKAGE = mfem
ABOUT = Lightweight, general, scalable C++ library for finite element methods
URL = https://github.com/mfem/mfem
PACKAGEVERSION = 4.8
## 4.9
## 4.5.2 has breaking changes for CHyPS.

# both i64
HYPREVERSION = 2.31.0
HYPRELOAD = ${HYPREVERSION}
PETSCVERSION = 3.23
PETSCLOAD = ${PETSCVERSION}

MODE = mpi
BASICMODULES = blaslapack adios2 hypre/${HYPRELOAD} metis petsc/${PETSCLOAD}
MODULES = ${BASICMODULES}

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

##
## petsc problem:
## https://github.com/mfem/mfem/issues/4526
##

WIDTH = 32
CMAKEFLAGS = -Wno-dev \
    -D CMAKE_BUILD_TYPE=Release \
    -D MFEM_USE_MPI=YES \
    -D MFEM_USE_PETSC=YES \
    -D PETSC_DIR=${PETSC_DIR} \
    -D PETSC_ARCH= \
    -D PETSC_EXECUTABLE_RUNS=${PETSC_DIR}/bin \
    \
    -D MFEM_USE_LAPACK=YES \
    -D LAPACK_LIBRARIES=${LAPACKLIB} \
    -D BLAS_LIBRARIES=${BLASLIB} \
    -D MFEM_DEBUG=NO \
    -D HYPRE_DIR=${TACC_HYPRE_DIR} \
    \
    -D IDXTYPEWIDTH=${WIDTH} -D REALTYPEWIDTH=${WIDTH} \
    -D CMAKE_C_FLAGS='-D IDXTYPEWIDTH=${WIDTH} -D REALTYPEWIDTH=${WIDTH}' \
    -D MFEM_USE_METIS_5=ON -D METIS_DIR=${TACC_METIS_DIR} \
    -D MFEM_USE_ADIOS2=ON -D ADIOS2_DIR=${TACC_ADIOS2_DIR}
nopetsc = \
    -D PETSC_LIBRARIES=${PETSC_DIR}/lib \
    -D PETSC_INCLUDES=${PETSC_DIR}/include \
    -D PETSC_EXECUTABLE_RUNS=${PETSC_DIR}/bin \
    -D DOXYGEN_EXECUTABLE=${TACC_DOXYGEN_BIN}/doxygen \
    -D foo=bar
CMAKEPREFIXPATHSET=1
CPTOINSTALLDIR = examples miniapps

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make cpu32 cpu64, gpu, (cpu is cpu32)"
	@echo "make default_install : cpu version"
.PHONY: default_install cpu gpu
default_install : cpu
cpu : cpu32
cpu32 : 
	@make --no-print-directory configure build public \
	    WIDTH=32
## https://github.com/mfem/mfem/issues/4844
cpu64 :
	@make --no-print-directory configure build public \
	    WIDTH=64 \
	    INSTALLVARIANT=i64 MODULEVERSIONEXTRA=i64 \
	    HYPRELOAD=${HYPREVERSION}-i64 PETSCLOAD=${PETSCVERSION}-i64

## https://github.com/mfem/mfem/issues/4527
CUDAFLAGS = -allow-unsupported-compiler -Xcompiler -std=gnu++14
gpu :
	@export CUDAFLAGS="${CUDAFLAGS}" \
	 && export NVCC_APPEND_FLAGS="${CUDAFLAGS}" \
	 && make --no-print-directory configure build public \
	    WIDTH=64 \
	    MODULES="${BASICMODULES} cuda" \
	    USEGPU=YES \
	    NVCCFLAGS="-D CMAKE_CUDA_FLAGS='${CUDAFLAGS}' " \
	    INSTALLVARIANT=cuda MODULEVERSIONEXTRA=cuda

##
## Module
##
DEPENDSONCURRENT = hypre

DOWNLOADURL = https://github.com/mfem/mfem/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = git@github.com:mfem/mfem.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean
