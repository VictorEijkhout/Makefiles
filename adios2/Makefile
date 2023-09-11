################################################################
####
#### Makefile for Adios2 installation
####
################################################################

PACKAGE = adios2
URL = https://github.com/adios2/Adios2
PACKAGEVERSION = git

MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

## PKGCONFIGLIB = pkgconfig

##
## cmake 
##
CMAKEFLAGS = -D CMAKE_BUILD_TYPE=Release
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

GITREPO = git@github.com:ornladios/ADIOS2.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

