################################################################
####
#### Makefile for Kokkos installation
####
################################################################

PACKAGE = KOKKOS
URL = https://github.com/kokkos/kokkos
PACKAGEVERSION = git
MODE = seq

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
    -D Kokkos_ENABLE_SERIAL=ON -D Kokkos_ENABLE_OPENMP=ON 
## VARSPROCESS = | sed -e 's/lib/lib64/'
## incompatible with openmp -D Kokkos_ENABLE_THREADS=ON
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
LIBDIR = lib64

info ::
	@echo "make cpu ( configure build )"
	@echo "make cuda ( configure_cuda build_cuda ) : cuda version in separate install dir"
.PHONY: cpu cuda configure_cuda build_cuda
cpu : configure build
# VOLTA70 ?
cuda : configure_cuda build_cuda
configure_cuda :
	make --no-print-directory configure \
	    CMAKEFLAGS="-D Kokkos_ENABLE_SERIAL=ON -D Kokkos_ENABLE_CUDA=ON -D Kokkos_ARCH_TURING75=ON -D Kokkos_ENABLE_CUDA_LAMBDA=ON" \
	    INSTALLEXT=cuda MODULES=cuda
build_cuda :
	make --no-print-directory build \
	    INSTALLEXT=cuda

GITREPO = git@github.com:kokkos/kokkos.git
include ${MAKEINCLUDES}/Make.git

