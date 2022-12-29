################################################################
####
#### Makefile for Kokkos installation
####
#### https://github.com/kokkos
####
################################################################

PACKAGE = KOKKOS
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = git
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODE = seq
MODULES = cuda

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
    -D Kokkos_ENABLE_SERIAL=ON -D Kokkos_ENABLE_OPENMP=ON 
## VARSPROCESS = | sed -e 's/lib/lib64/'
## incompatible with openmp -D Kokkos_ENABLE_THREADS=ON
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
LIBDIR = lib64

info ::
	@echo "make configure_cuda build_cuda : cuda version in separate install dir"
.PHONY: configure_cuda build_cuda
# VOLTA70 ?
configure_cuda :
	make --no-print-directory configure \
	    CMAKEFLAGS="-D Kokkos_ENABLE_SERIAL=ON -D Kokkos_ENABLE_CUDA=ON -D Kokkos_ARCH_TURING75=ON -D Kokkos_ENABLE_CUDA_LAMBDA=ON" \
	    INSTALLEXT=cuda
build_cuda :
	make --no-print-directory build \
	    INSTALLEXT=cuda

GITREPO = git@github.com:kokkos/kokkos.git
include ${MAKEINCLUDES}/Make.git

