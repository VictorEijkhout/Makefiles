################################################################
####
#### Makefile for Metis as patched by scivision installation
####
#### https://github.com/scivision/METIS
#### this the 5.1.0.3 version only
#### there are big problems with the Karypis 5.2 version
####
################################################################

PACKAGE = metis
PACKAGEVERSION = 5.1.0.3
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/scivision/METIS/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

