################################################################
####
#### Makefile for P4est installation
####
################################################################

PACKAGE = p4est

# PACKAGEVERSION = 2.8.5
# TGZURL = https://p4est.github.io/release/p4est-${PACKAGEVERSION}.tar.gz

PACKAGEVERSION = 2.8.5.5-9ddbb
TGZURL = https://p4est.github.io/release/p4est-${PACKAGEVERSION}.tar.gz

MODE = mpi
MODULES = petsc

include ${MAKEINCLUDES}/Make.info

CONFIGUREFLAGS=--enable-mpi
## CMAKEFLAGS=--enable-mpi
PKGCONFIG = lib/pkgconfig

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

##TGZURL = https://github.com/p4est/p4est.github.io/blob/master/release/p4est-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/cburstedde/p4est.git
BRANCH = master
SUBMODULE = yes
PULLAFTERACTIONS = ./bootstrap
include ${MAKEINCLUDES}/Make.git
