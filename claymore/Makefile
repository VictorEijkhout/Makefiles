################################################################
####
#### Makefile for Claymore installation
####
################################################################

PACKAGE = claymore
URL = https://github.com/Javiermro/claymore_JLM
PACKAGEVERSION = git

MODE = mpi cuda
MODULES = 

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
## https://github.com/claymore/claymore/issues/4526
##

CMAKEFLAGS = 
CMAKEPREFIXPATHSET=1

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make cpu, gpu"
	@echo "make default_install : cpu version"
.PHONY: default_install gpu
default_install : gpu
CUDAFLAGS = -allow-unsupported-compiler -Xcompiler -std=gnu++14
gpu :
	@export CUDAFLAGS="${CUDAFLAGS}" \
	 && export NVCC_APPEND_FLAGS="${CUDAFLAGS}" \
	 && make --no-print-directory configure build \
	    USEGPU=YES \
	    NVCCFLAGS="-D CMAKE_CUDA_FLAGS='${CUDAFLAGS}' "

GITREPO = https://github.com/Javiermro/claymore_JLM.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean
