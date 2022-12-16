################################################################
####
#### Makefile for P4est installation
####
################################################################

PACKAGE = p4est
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 2.8
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODE = mpi
MODULES = petsc

include ${MAKEINCLUDES}/Make.info

CONFIGUREFLAGS="LIBS=-mkl --enable-mpi"
PKGCONFIGSET = lib/pkgconfig
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install

TGZURL = https://p4est.github.io/release/p4est-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/cburstedde/p4est.git
SUBMODULE = yes
PULLAFTERACTIONS = ./bootstrap
include ${MAKEINCLUDES}/Make.git
