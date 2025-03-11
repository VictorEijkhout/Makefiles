################################################################
####
#### Makefile for Openmm installation
####
################################################################

PACKAGE = openmmcuda
URL = https://github.com/openmm/openmm
PACKAGEVERSION = 8.2.0

MODE = mpi
BASICMODULES = 
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
    -D OPENMM_USE_CUDA=${USEGPU} \
    -D CMAKE_CUDA_FLAGS="-std=c++14" \
    -D LAPACK_LIBRARIES=${LAPACKLIB} \
    -D BLAS_LIBRARIES=${BLASLIB} 
## CMAKEPREFIXPATHSET=1

include ${MAKEINCLUDES}/Make.cmake
info ::
	@echo "make config-help"
config-help :
	@cd ${STOCKYARD}/openmm/openmm-${PACKAGEVERSION} \
	 && open docs-source/usersguide/application/01_getting_started.rst
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : gpu version"
.PHONY: default_install gpu
default_install gpu :
	@export CUDAFLAGS="${CUDAFLAGS}" \
	 && export NVCC_APPEND_FLAGS="${CUDAFLAGS}" \
	 && make --no-print-directory configure build 

TGZURL = https://github.com/openmm/openmm/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = git@github.com:openmm/openmm.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean
