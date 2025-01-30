################################################################
####
#### Makefile for Entity installation
####

################################################################

PACKAGE = entity
URL = https://github.com/entity-toolkit/entity
PACKAGEVERSION = 1.1.1

MODE = mpi
MODULES_BASIC = phdf5
MODULES = ${MODULES_BASIC}

include ${MAKEINCLUDES}/Make.info

##
## cmake version
##
CMAKEFLAGS_BASIC =  \
-D output=ON -D mpi=ON
CMAKEFLAGS = ${CMAKEFLAGS_BASIC}

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : cpu, also gpu"
.PHONY: default_install cpu gpu
default_install : cpu
cpu : configure build
gpu :
	@make --no-print-directory configure build \
	  MODULES="${MODULES_BASIC} cuda" \
	  CMAKEFLAGS="${CMAKEFLAGS_BASIC} -D Kokkos_ENABLE_CUDA=ON -D Kokkos_ARCH_HOPPER90=ON"

TGZURL = https://github.com/entity-toolkit/entity/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = 
include ${MAKEINCLUDES}/Make.git

