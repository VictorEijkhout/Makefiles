################################################################
####
#### Makefile for Kokkos installation
####
################################################################

PACKAGE = KOKKOS
URL = https://github.com/kokkos/kokkos \
    https://kokkos.github.io/kokkos-core-wiki/building.html
PACKAGEVERSION = 4.6.01
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
	@echo "make cpu ( configure build )"
	@echo "make gpu ( configure_cuda build_cuda ) : cuda version in separate install dir"
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
.PHONY: gpu configure_cuda 
gpu : 
	@make --no-print-directory \
	    configure_cuda build \
	    MODULES=cuda DEPENDSON=cuda/12 \
	    INSTALLVARIANT=cuda INSTALLEXT=cuda
ziplog_info :
	@echo "Please explicit ziplog-cpu or ziplog-gpu:"
	@echo make ziplog INSTALLVARIANT=omp INSTALLEXT=omp
	@echo make ziplog INSTALLVARIANT=cuda INSTALLEXT=cuda
## https://www.gnu.org/software/make/manual/html_node/Overriding-Makefiles.html

# VOLTA70 ?
# https://kokkos.org/kokkos-core-wiki/get-started/configuration-guide.html
configure_cuda :
	@if [ -z "${TACC_CXX}" ] ; then \
	    echo "No TACC_CXX defined" && exit 1 ; fi 
	@export NVCC_WRAPPER_DEFAULT_COMPILER=${TACC_CXX} \
	 && case ${TACC_SYSTEM} in \
	    ( vista) GPUTYPE=HOPPER90 ;; ( * ) GPUTYPE=TURING75 ;; esac \
	 && make --no-print-directory configure \
	    CMAKEFLAGS=" \
	      -D Kokkos_ENABLE_SERIAL=ON \
	      -D Kokkos_ENABLE_CUDA=ON -D Kokkos_ENABLE_CUDA_LAMBDA=ON \
	      -D Kokkos_ARCH_$${GPUTYPE}=ON \
	      " \
	    MODULES=cuda  INSTALLVARIANT=cuda
## NO LONGER (at least on vista)
##	      -D CMAKE_CXX_COMPILER=${SRCPATH}/bin/nvcc_wrapper \
## 	      -D KOKKOS_CUDA_OPTIONS=-allow-unsupported-compiler \
##      
##   -D CMAKE_CXX_FLAGS=--gcc-toolchain=/opt/apps/gcc/12.2.0

.PHONY: sycl configure_sycl 
sycl : 
	@make --no-print-directory PACKAGE_VERSION=${PACKAGE_VERSION} \
	    configure_sycl build public \
	    INSTALLVARIANT=sycl INSTALLEXT=sycl
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
	    INSTALLVARIANT=sycl

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
