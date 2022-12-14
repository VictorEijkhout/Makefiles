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
include ${MAKEINCLUDES}/Make.build

TGZURL = https://github.com/cburstedde/p4est/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
