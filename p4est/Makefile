################################################################
####
#### Makefile for P4est installation
####
################################################################

PACKAGE = p4est
ABOUT = Parallel AMR on Forests of Octrees
URL = http://p4est.org/
CODEURL = https://github.com/cburstedde/p4est
PACKAGEVERSION = 2.8.7

MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

# hack'y solution for:
# https://github.com/cburstedde/p4est/issues/351
PKGCONFIG = lib/pkgconfig
CMAKEPREFIXPATHSET = 1

#ifeq "${TACC_FAMILY_COMPILER}" "gcc"
ifeq "foo" "foo"

##
## cmake version
##
CMAKEFLAGS = \
    -D SC_ENABLE_MPI=ON \
    -D P4EST_ENABLE_MPI=ON \
    -D P4EST_BUILD_EXAMPLES=OFF
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

else 

##
## configure version
##
CONFIGUREFLAGS=--enable-mpi
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install

endif

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

##
## download
##
DOWNLOADURL = https://github.com/cburstedde/p4est/releases/download/v${PACKAGEVERSION}/p4est-${PACKAGEVERSION}.tar.gz
## https://p4est.github.io/release/p4est-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/cburstedde/p4est.git
BRANCH = master
SUBMODULE = yes
PULLAFTERACTIONS = ./bootstrap
include ${MAKEINCLUDES}/Make.git
