################################################################
####
#### Makefile for Kokkos installation
####
################################################################

PACKAGE = KOKKOS
URL = https://github.com/kokkos/kokkos \
    https://kokkos.github.io/kokkos-core-wiki/building.html
PACKAGEVERSION = 4.1.00
MODE = seq

include ${MAKEINCLUDES}/Make.info

CMAKE_BASIC_FLAGS = \
    -D CMAKE_VERBOSE_MAKEFILE=ON 

## VARSPROCESS = | sed -e 's/lib/lib64/'
## incompatible with openmp -D Kokkos_ENABLE_THREADS=ON
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
LIBDIR = lib64

info ::
	@echo "make cpu ( configure build )"
	@echo "make cuda ( configure_cuda build_cuda ) : cuda version in separate install dir"
.PHONY: cpu gpu configure_cuda build_cuda
CMAKEFLAGS = \
    ${CMAKE_BASIC_FLAGS} \
    -D Kokkos_ENABLE_SERIAL=ON -D Kokkos_ENABLE_OPENMP=ON 
cpu : configure build
# VOLTA70 ?
gpu : configure_cuda build_cuda
configure_cuda :
	@if [ -z "${CXX}" ] ; then \
	    echo "No CXX defined" && exit 1 ; fi 
	@export NVCC_WRAPPER_DEFAULT_COMPILER=${CXX} \
	 && make --no-print-directory configure \
	    CMAKEFLAGS="${CMAKE_BASIC_FLAGS} \
	      -D CMAKE_CXX_COMPILER=${PACKAGEROOT}/kokkos/kokkos-${PACKAGEVERSION}/bin/nvcc_wrapper \
	      -D CMAKE_CXX_FLAGS=-allow-unsupported-compiler \
	      -D Kokkos_ENABLE_SERIAL=ON -D Kokkos_ENABLE_CUDA=ON -D Kokkos_ARCH_TURING75=ON -D Kokkos_ENABLE_CUDA_LAMBDA=ON \
	      -D KOKKOS_CUDA_OPTIONS=-allow-unsupported-compiler \
	      " \
	    INSTALLEXT=cuda MODULES=cuda
build_cuda :
	make --no-print-directory build \
	    INSTALLEXT=cuda

TGZURL = https://github.com/kokkos/kokkos/archive/refs/tags/4.1.00.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/kokkos/kokkos.git
include ${MAKEINCLUDES}/Make.git

