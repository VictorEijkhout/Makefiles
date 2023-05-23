################################################################
####
#### Makefile for Netcdf installation
####
################################################################

PACKAGE = SUITESPARSE
PACKAGEVERSION = 7.0.1
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

