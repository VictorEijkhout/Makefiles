################################################################
####
#### Makefile for Enzo installation
####
################################################################

PACKAGE = enzo
URL = http://enzo.org/ https://github.com/cburstedde/enzo
PACKAGEVERSION = git

MODE = mpi
MODULES = boost charmpp phdf5

include ${MAKEINCLUDES}/Make.info

##
## cmake version
##
CMAKEFLAGS = \
    -D CHARM_ROOT=${TACC_CHARMPP_DIR} \
    -D USE_GRACKLE=OFF \
    -D USE_DOUBLE_PREC=OFF \
    \
    -D use_jemalloc=OFF \
    -D balance=OFF  \
    -D darshan=OFF \
    -D memory=OFF \
    -D smp=OFF \
    -D use_projections=OFF
HASBIN = 1

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

# TGZURL = 
# include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/enzo-project/enzo-e.git
## BRANCH = master
## SUBMODULE = yes
## PULLAFTERACTIONS = ./bootstrap
include ${MAKEINCLUDES}/Make.git
