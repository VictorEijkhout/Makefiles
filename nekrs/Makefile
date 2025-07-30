################################################################
####
#### Makefile for NekRS installation
####
################################################################

PACKAGE = nekRS
URL = https://github.com/Nek5000/nekRS
PACKAGEVERSION = 23.0

MODE = mpi
MODULES_BASIC = 
MODULES = ${MODULES_BASIC} 

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
CMAKEFLAGS_BASIC = \
    -D USE_OPTIMIZED_BLAS=ON \
    -D BLAS_LIBRARIES='${BLASLIB}' \
    -D USE_OPTIMIZED_LAPACK=ON \
    -D LAPACK_LIBRARIES='${LAPACKLIB}'
CMAKEFLAGS = ${CMAKEFLAGS_BASIC}
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
	  MODULES="${MODULES_BASIC} cuda" INSTALLVARIANT=cuda MODULEVERSIONEXTRA=cuda \
	  CMAKEFLAGS="${CMAKEFLAGS_BASIC} -D NEKRS_USE_CUDA:BOOL=ON"

##
## Download
##
TGZURL = https://github.com/Nek5000/nekRS/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

##
## Clean
##
include ${MAKEINCLUDES}/Make.clean

