################################################################
####
#### Makefile for P4est installation
####
################################################################

PACKAGE = p4est
URL = http://p4est.org/ https://github.com/cburstedde/p4est
PACKAGEVERSION = 2.8.6
## PACKAGEVERSION = 2.8.5.5-9ddbb

MODE = mpi
MODULES = petsc

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig

##
## configure version
##
CONFIGUREFLAGS=--enable-mpi
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install

##
## cmake version
##
# CMAKEFLAGS=--enable-mpi
# include ${MAKEINCLUDES}/Make.cmake
# include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://p4est.github.io/release/p4est-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/cburstedde/p4est.git
BRANCH = master
SUBMODULE = yes
PULLAFTERACTIONS = ./bootstrap
include ${MAKEINCLUDES}/Make.git
