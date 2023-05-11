################################################################
####
#### Makefile for Mdspan installation
####
#### https://github.com/kokkos/mdspan
####
################################################################

PACKAGE = MDSPAN
PACKAGEVERSION = git
MODE = seq

include ${MAKEINCLUDES}/Make.info

## CMAKEFLAGS = 
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
LIBDIR = lib64

GITREPO = git@github.com:kokkos/mdspan.git
include ${MAKEINCLUDES}/Make.git

