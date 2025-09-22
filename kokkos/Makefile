################################################################
####
#### Makefile for Kokkos installation
####
################################################################

PACKAGE = KOKKOS
ABOUT = C++ Performance Portability Programming Ecosystem
URL = https://github.com/kokkos/kokkos
DOCURL = https://kokkos.github.io/kokkos-core-wiki/building.html
PACKAGEVERSION = 4.6.01
# 4.7.01 is out
MODE = seq
MODULES = 
MODULE_MINIMUM_VERSION_cuda = 12

include ${MAKEINCLUDES}/Make.info

CMAKE_BASIC_FLAGS = 

## VARSPROCESS = | sed -e 's/lib/lib64/'
## incompatible with openmp -D Kokkos_ENABLE_THREADS=ON
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
DEPENDSON = 

info ::
	@echo "make cpu  ( configure build )"
	@echo "make cuda ( configure_cuda build_cuda ) : cuda version in separate install dir"
	@echo "make sycl ( configure_sycl build_sycl ) : sycl version in separate install dir"
.PHONY: cpu 
.PHONY: default_install
info ::
	@echo "make default_install : omp only, do make gpu for cuda backend"
default_install : cpu
CMAKEFLAGS = \
    -D Kokkos_ENABLE_SERIAL=ON -D Kokkos_ENABLE_OPENMP=ON 
HASBIN = 1
## we need this so that we can override it in jail
CMAKEPREFIXPATHSET = 1
SRCPATH = ${PACKAGEROOT}/kokkos/kokkos-${PACKAGEVERSION}
cpu :
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

.PHONY: cuda gpu configure_cuda 
gpu cuda : configure_cuda
	@make --no-print-directory \
	    build public \
	    MODULES=cuda DEPENDSON=cuda/12 \
	    INSTALLVARIANT=cuda INSTALLEXT=cuda \
	    EXTRAVARS="Kokkos_ENABLE_CUDA=ON Kokkos_ARCH_${COMPUTE}=ON Kokkos_ARCH_${HOST}=ON"

# VOLTA70 ?
# https://kokkos.org/kokkos-core-wiki/get-started/configuration-guide.html
configure_cuda :
	@if [ -z "${TACC_CXX}" ] ; then \
	    echo "No TACC_CXX defined" && exit 1 ; fi 
	@export NVCC_WRAPPER_DEFAULT_COMPILER=${TACC_CXX} \
	 && make --no-print-directory configure \
	    CMAKEFLAGS=" \
	      -D Kokkos_ENABLE_SERIAL=ON \
	      -D Kokkos_ENABLE_CUDA=ON -D Kokkos_ENABLE_CUDA_LAMBDA=ON \
	      -D Kokkos_ARCH_${COMPUTE}=ON \
	      " \
	    INSTALLVARIANT=cuda INSTALLEXT=cuda \
	    MODULES=cuda

##
## SYCL
##
.PHONY: sycl configure_sycl 
sycl : 
	@make --no-print-directory PACKAGE_VERSION=${PACKAGE_VERSION} \
	    configure_sycl build public \
	    INSTALLVARIANT=sycl INSTALLEXT=sycl \
	    EXTRAVARS="Kokkos_ENABLE_SYCL=ON Kokkos_ARCH_INTEL_PVC=ON Kokkos_ARCH_${HOST}=ON"

# VOLTA70 ?
configure_sycl :
	@if [ -z "${TACC_CXX}" ] ; then \
	    echo "No TACC_CXX defined" && exit 1 ; fi 
	@export NVCC_WRAPPER_DEFAULT_COMPILER=${TACC_CXX} \
	 && make --no-print-directory configure \
	    CMAKEFLAGS=" \
	      -D CMAKE_CXX_COMPILER=icpx \
	      -D Kokkos_ENABLE_SERIAL=ON -D Kokkos_ENABLE_SYCL=ON \
	      " \
	    INSTALLVARIANT=sycl INSTALLEXT=sycl

info ::
	@echo "make kokkoscore INSTALLVARIANT=..."
.PHONY: kokkoscore
kokkoscore :
	@source ${MAKEINCLUDES}/names.sh \
	 && export INSTALLROOT=${INSTALLROOT} \
	 && setnames \
	    "${PACKAGE}"  "${PACKAGEVERSION}"  "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}"  "${INSTALLROOT}"  "${INSTALLEXT}"  "$${petscarch}" \
	    "${HOMEDIR}"  "${BUILDDIRROOT}"  "${MODE}" \
	    "${PREFIXOPTION}"  "${PREFIXEXTRA}" \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 \
	 && echo $${builddir} && ls $${builddir} \
	 && cat $${builddir}/KokkosCore_config.h

TGZURL = https://github.com/kokkos/kokkos/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/kokkos/kokkos.git
include ${MAKEINCLUDES}/Make.git

clean ::
	@rm -f cuda_install.o*
