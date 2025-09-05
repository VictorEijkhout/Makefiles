################################################################
####
#### Makefile for CANTERA installation
####
################################################################

PACKAGE = CANTERA
ABOUT = 
PACKAGEVERSION = 3.1.0
URL = https://github.com/Cantera/cantera
MODE = seq
MODULES = boost mkl

include ${MAKEINCLUDES}/Make.info

##
## CMake
##

CMAKESUBDIR = interfaces/python_sdist

# Blas / Lapack
ifeq "${TACC_SYSTEM}" "vista"
  BLASLIB=${TACC_NVPL_LIB}/libnvpl_blas_lp64_seq.so 
  LAPACKLIB=${TACC_NVPL_LIB}/libnvpl_lapack_lp64_seq.so
else
  # mkl_intel_lp64;mkl_sequential;mkl_core;iomp5;pthread"
  BLASLIB=${TACC_MKL_LIB}/libmkl_sequential.so ${TACC_MKL_LIB}/libmkl_core.so 
  LAPACKLIB=${TACC_MKL_LIB}/libmkl_intel_lp64.so ${TACC_MKL_LIB}/libmkl_sequential.so 
endif

CMAKEFLAGS = -Wno-dev \
    -D LAPACK_LIBRARIES='${LAPACKLIB}' \
    -D BLAS_LIBRARIES='${BLASLIB}'

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build public

####
#### Download
####
TGZURL = https://github.com/Cantera/cantera/releases/download/v${PACKAGEVERSION}/cantera-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/Cantera/cantera.git
include ${MAKEINCLUDES}/Make.git
