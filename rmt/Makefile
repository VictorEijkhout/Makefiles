################################################################
####
#### Makefile for Rmt installation
####
################################################################

PACKAGE = rmt
URL = https://gitlab.com/Uk-amor/RMT/rmt
PACKAGEVERSION = git

MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

## PKGCONFIGLIB = pkgconfig

##
## cmake version
##
# CMAKEFLAGS=--enable-mpi
include ${MAKEINCLUDES}/Make.cmake
CMAKESUBDIR = source
include ${MAKEINCLUDES}/Make.cbuild

GITREPO = https://gitlab.com/Uk-amor/RMT/rmt.git
include ${MAKEINCLUDES}/Make.git
