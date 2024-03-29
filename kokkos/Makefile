################################################################
####
#### Makefile for Kokkos installation
####
################################################################

PACKAGE = KOKKOS
URL = https://github.com/kokkos/kokkos \
    https://kokkos.github.io/kokkos-core-wiki/building.html
PACKAGEVERSION = 4.2.00
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
.PHONY: cpu gpu configure_cuda
.PHONY: default_install
info ::
	@echo "make default_install : omp only, do make gpu for cuda backend"
default_install : cpu
CMAKEFLAGS = \
    -D Kokkos_ENABLE_SERIAL=ON -D Kokkos_ENABLE_OPENMP=ON 
## we need this so that we can override it in jail
CMAKEPREFIXPATHSET = 1
SRCPATH = ${PACKAGEROOT}/kokkos/kokkos-${PACKAGEVERSION}
cpu :
	@make --no-print-directory PACKAGE_VERSION=${PACKAGE_VERSION} \
	    configure build public \
	    INSTALLVARIANT=omp INSTALLEXT=omp
gpu : 
	@make --no-print-directory PACKAGE_VERSION=${PACKAGE_VERSION} \
	    configure_cuda build public \
	    MODULES=cuda DEPENDSON=cuda/12 \
	    VARIANT=cuda INSTALLEXT=cuda
# VOLTA70 ?
configure_cuda :
	@if [ -z "${TACC_CXX}" ] ; then \
	    echo "No TACC_CXX defined" && exit 1 ; fi 
	@export NVCC_WRAPPER_DEFAULT_COMPILER=${TACC_CXX} \
	 && make --no-print-directory configure \
	    CMAKEFLAGS=" \
	      -D CMAKE_CXX_COMPILER=${SRCPATH}/bin/nvcc_wrapper \
	      -D Kokkos_ENABLE_SERIAL=ON -D Kokkos_ENABLE_CUDA=ON \
	      -D Kokkos_ARCH_TURING75=ON -D Kokkos_ENABLE_CUDA_LAMBDA=ON \
	      -D KOKKOS_CUDA_OPTIONS=-allow-unsupported-compiler \
	      " \
	    MODULES=cuda  VARIANT=cuda
## CMAKE_CXX_FLAGS=-allow-unsupported-compiler
##   -D CMAKE_CXX_FLAGS=--gcc-toolchain=/opt/apps/gcc/12.2.0

info ::
	@echo "make kokkoscore INSTALLVARIANT=..."
.PHONY: kokkoscore
kokkoscore :
	@source ${MAKEINCLUDES}/names.sh \
	 && export INSTALLROOT=${INSTALLROOT} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 \
	 && echo $${builddir} && ls $${builddir} \
	 && cat $${builddir}/KokkosCore_config.h

TGZURL = https://github.com/kokkos/kokkos/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/kokkos/kokkos.git
include ${MAKEINCLUDES}/Make.git

