################################################################
####
#### Makefile for Entity installation
####
################################################################

PACKAGE = entity
ABOUT = New generation astrophysical plasma simulation code with CPU/GPU portability
URL = https://github.com/entity-toolkit/entity
PACKAGEVERSION = 1.2.1
# 1.2.5

MODE = mpi
MODULES_BASIC = kokkos adios2 phdf5
MODULES = ${MODULES_BASIC}

include ${MAKEINCLUDES}/Make.info

##
## cmake version
##
# see `pgens' directory for choices
CMAKEFLAGS_BASIC =  \
    -D output=ON -D mpi=ON \
    -D pgen=shock \
    -D HDF5_ROOT=${TACC_PHDF5_DIR}
CMAKEFLAGS = ${CMAKEFLAGS_BASIC}

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : cpu, also gpu"
.PHONY: default_install cpu gpu
default_install : cpu
HASBIN = 1
NOINC = 1
NOLIB = 1
cpu :
	@export HDF5_ROOT=${TACC_PHDF5_DIR} \
	 && ${MAKE} --no-print-directory configure build 
gpu :
	@export HDF5_ROOT=${TACC_PHDF5_DIR} \
	 && export NVCC_APPEND_FLAGS='-allow-unsupported-compiler' \
	 && make --no-print-directory configure build \
	  MODULES="${MODULES_BASIC} cuda" INSTALLVARIANT=cuda MODULEVERSIONEXTRA=cuda 
	  CMAKEFLAGS="${CMAKEFLAGS_BASIC} -D Kokkos_ENABLE_CUDA=ON "
## -D Kokkos_ARCH_HOPPER90=ON"

TGZURL = https://github.com/entity-toolkit/entity/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/entity-toolkit/entity.git
include ${MAKEINCLUDES}/Make.git
