################################################################
####
#### Makefile for NekRS installation
####
################################################################

PACKAGE = nekRS
ABOUT = Next generation fast and scalable CFD code
URL = https://github.com/Nek5000/nekRS
PACKAGEVERSION = 26.0

MODE = mpi
MODULES_BASIC = 
MODULES = ${MODULES_BASIC} 
MODULES_gpu = "${MODULES_BASIC} cuda"
include ${MAKEINCLUDES}/Make.info

##
## Cmake 
##
ifeq "${TACC_SYSTEM}" "vista"
  BLASLIB=${TACC_NVPL_LIB}/libnvpl_blas_lp64_seq.so 
  LAPACKLIB=${TACC_NVPL_LIB}/libnvpl_lapack_lp64_seq.so
else
  # mkl_intel_lp64;mkl_sequential;mkl_core;iomp5;pthread"
  BLASLIB=${TACC_MKL_LIB}/libmkl_intel_lp64.so ${TACC_MKL_LIB}/libmkl_core.so 
  LAPACKLIB=${BLASLIB}
endif
OPTIMIZEDBLAS = OFF

USECUDA = OFF
# cuda fails on hypre:
# https://github.com/Nek5000/nekRS/issues/648
CMAKEFLAGS = \
    -D USE_OPTIMIZED_BLAS=${OPTIMIZEDBLAS} \
    -D BLAS_LIBRARIES='${BLASLIB}' \
    -D USE_OPTIMIZED_LAPACK=${OPTIMIZEDBLAS} \
    -D LAPACK_LIBRARIES='${LAPACKLIB}' \
    -D NEKRS_USE_CUDA:BOOL=${USECUDA} \
    -D OCCA_ENABLE_CUDA=${USECUDA} \
    -D NEKRS_GPU_MPI=OFF
HASBIN = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : cpu; also: gpu"
.PHONY: default_install cpu gpu
default_install : cpu
cpu : configure build
gpu :
	@make --no-print-directory configure build \
	  MODULES="${MODULES_gpu}" INSTALLVARIANT=cuda MODULEVERSIONEXTRA=cuda 

##
## Download
##
DOWNLOADURL = https://github.com/Nek5000/nekRS/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
BRANCH = next
# BOOTSTRAP = find . -name polyfill.hpp -exec sed -i -e 's/ uint/ std::uint/g' {} \;
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/Nek5000/nekRS.git
include ${MAKEINCLUDES}/Make.git

##
## Clean
##
include ${MAKEINCLUDES}/Make.clean

