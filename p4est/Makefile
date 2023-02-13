################################################################
####
#### Makefile for P4est installation
####
################################################################

PACKAGE = p4est
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 2.8
MODE = mpi
MODULES = petsc

include ${MAKEINCLUDES}/Make.info

CONFIGUREFLAGS=--enable-mpi
## CMAKEFLAGS=--enable-mpi
PKGCONFIGSET = lib/pkgconfig

##
## configure version
##
#LIBS = -lirc
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
##
## cmake version
##
# include ${MAKEINCLUDES}/Make.cmake
# include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://p4est.github.io/release/p4est-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/cburstedde/p4est.git
SUBMODULE = yes
PULLAFTERACTIONS = ./bootstrap
include ${MAKEINCLUDES}/Make.git
