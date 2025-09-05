################################################################
####
#### Makefile for Ginkgo installation
####
################################################################

PACKAGE = GINKGO
ABOUT = High-Performance Linear Algebra Library for Manycore Systems
URL = https://github.com/ginkgo-project/ginkgo
DOCURL = https://github.com/ginkgo-project/ginkgo/blob/develop/INSTALL.md
PACKAGEVERSION = 1.10.0
MODE = mpi
MODULES = 
MODULE_MINIMUM_VERSION_cuda = 12

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
DEPENDSON = 

info ::
	@echo "make cpu  ( configure build )"
	@echo "make cuda ( configure_cuda build_cuda ) : cuda version in separate install dir"
	@echo "make sycl ( configure_sycl build_sycl ) : sycl version in separate install dir"
info ::
	@echo "make default_install, omp, cpu : omp only, do make gpu for cuda backend"
.PHONY: default_install
default_install : cpu
# disabling anything involving linking
# https://github.com/ginkgo-project/ginkgo/issues/1928
LINKEDSTUFF = ON
CMAKEFLAGS = \
    -D GINKGO_BUILD_BENCHMARKS=${LINKEDSTUFF} \
    -D GINKGO_BUILD_EXAMPLES=${LINKEDSTUFF} \
    -D GINKGO_BUILD_TESTS=${LINKEDSTUFF} \
    -D GINKGO_BUILD_OMP=ON
HASBIN = 1
## we need this so that we can override it in jail
CMAKEPREFIXPATHSET = 1
SRCPATH = ${PACKAGEROOT}/ginkgo/ginkgo-${PACKAGEVERSION}
.PHONY: cpu omp
cpu omp :
	@make --no-print-directory \
	    configure build \
	    INSTALLVARIANT=omp INSTALLEXT=omp

ziplog_info :
	@echo "Please explicit ziplog-cpu or ziplog-gpu:"
	@echo make ziplog INSTALLVARIANT=omp INSTALLEXT=omp
	@echo make ziplog INSTALLVARIANT=cuda INSTALLEXT=cuda
## https://www.gnu.org/software/make/manual/html_node/Overriding-Makefiles.html

HOST_fontera = CLX
HOST_ls6 = MILAN
HOST_stampede3 = SPR
HOST  = ${HOST_${TACC_SYSTEM}}

COMPUTE_frontera = TURING75
COMPUTE_ls6 = TURING75
COMPUTE_stampede3 = HOPPER90
COMPUTE = ${COMPUTE_${TACC_SYSTEM}}

##
## CUDA
##

.PHONY: cuda gpu configure_cuda build_cuda
gpu cuda : configure_cuda build_cuda
configure_cuda :
	@make --no-print-directory configure \
	    CMAKEFLAGS=" \
	        ${CMAKEFLAGS} \
	        -D GINKGO_BUILD_CUDA=ON \
	      " \
	    INSTALLVARIANT=cuda INSTALLEXT=cuda \
	    MODULES=cuda
build_cuda :
	@make --no-print-directory \
	    build public \
	    CMAKEFLAGS=" \
	        ${CMAKEFLAGS} \
	        -D GINKGO_BUILD_CUDA=ON \
	      " \
	    MODULES=cuda DEPENDSON=cuda/12 \
	    INSTALLVARIANT=cuda INSTALLEXT=cuda

##
## SYCL
##
.PHONY: sycl configure_sycl 
sycl : 
	@make --no-print-directory PACKAGE_VERSION=${PACKAGE_VERSION} \
	    configure_sycl build public \
	    INSTALLVARIANT=sycl INSTALLEXT=sycl \
	    EXTRAVARS="Ginkgo_ENABLE_SYCL=ON Ginkgo_ARCH_INTEL_PVC=ON Ginkgo_ARCH_${HOST}=ON"

# VOLTA70 ?
configure_sycl :
	@if [ -z "${TACC_CXX}" ] ; then \
	    echo "No TACC_CXX defined" && exit 1 ; fi 
	@export NVCC_WRAPPER_DEFAULT_COMPILER=${TACC_CXX} \
	 && make --no-print-directory configure \
	    CMAKEFLAGS=" \
	      -D CMAKE_CXX_COMPILER=icpx \
	      -D Ginkgo_ENABLE_SERIAL=ON -D Ginkgo_ENABLE_SYCL=ON \
	      " \
	    INSTALLVARIANT=sycl INSTALLEXT=sycl

TGZURL = https://github.com/ginkgo-project/ginkgo/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

clean ::
	@rm -f cuda_install.o*
