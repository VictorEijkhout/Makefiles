################################################################
####
#### Makefile for Peridigm installation
####
#### needs a modern enough C++: gcc12 or so.
####
################################################################

PACKAGE = PERIDIGM
ABOUT = Peridigm peridynamics code.
URL = https://github.com/peridigm/peridigm
PACKAGEVERSION = git
MODE = mpi
MODULES = trilinos boost pnetcdf
# trilinos needs to be /13.0.0 ???

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
    -D CMAKE_BUILD_TYPE:STRING=Release \
    -D Trilinos_DIR:PATH=${TACC_TRILINOS_DIR} \
    -D Trilinos_INCLUDE_DIRS=${TACC_TRILINOS_INC}
ifeq "${TACC_FAMILY_COMPILER}" "gcc"
CMAKEFLAGS += \
    -D CMAKE_CXX_FLAGS=\"-fopenmp -std=c++17\"
else
CMAKEFLAGS += \
    -D CMAKE_CXX_FLAGS=\"-qopenmp -std=c++17\"
endif

#     -D Trilinos_INCLUDE_DIRS='${TACC_TRILINOS_INC};${TACC_NETCDF_INC}'

include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
##
## total hack because peridigm is not up with trilinos >2022
##
# CPTOSOURCE = src,${TACC_NETCDF_INC}/\*.h
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build

GITREPO = https://github.com/peridigm/peridigm.git
include ${MAKEINCLUDES}/Make.git
include ${MAKEINCLUDES}/Make.clean

