################################################################
####
#### Makefile for Enzo installation
####
################################################################

PACKAGE = enzo
URL = https://github.com/enzo-project/enzo-e
PACKAGEVERSION = git
MODE = mpi
MODULES = boost charmpp grackle hdf5

include ${MAKEINCLUDES}/Make.info

##
## cmake version
##
CMAKEFLAGS = \
    -D Enzo-E_CONFIG=linux_nvc \
    -D use_projections=OFF \
    -D USE_YT_BASED_TESTS=OFF \
    -D USE_DOUBLE_PREC=OFF \
    \
    -D CHARM_ROOT=${TACC_CHARMPP_DIR} \
    -D USE_GRACKLE=OFF \
    -D Grackle_ROOT=${TACC_GRACKLE_DIR} \
    -D Grackle_LIBRARY=${TACC_GRACKLE_LIB} \
    -D GRACKLE_INPUT_DATA_DIR=${TACC_GRACKLE_DIR}/grackle_data_files/input \
    \
    -D use_jemalloc=OFF \
    -D balance=OFF  \
    -D darshan=OFF \
    -D memory=OFF \
    -D smp=OFF
CPPSTANDARD = 20
BUILDSTATICLIBS	= 1
HASBIN = 1

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

# TGZURL = 
# include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/jobordner/enzo-e.git
BRANCH = port/tacc-grace
## SUBMODULE = yes
## PULLAFTERACTIONS = ./bootstrap
include ${MAKEINCLUDES}/Make.git
