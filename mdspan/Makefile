################################################################
####
#### Makefile for Mdspan installation
####
################################################################

PACKAGE = MDSPAN
URL = https://github.com/kokkos/mdspan
PACKAGEVERSION = git
MODE = seq

include ${MAKEINCLUDES}/Make.info

## CMAKEFLAGS = 
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
LIBDIR = lib64
.PHONY: default_install
default_install : configure build

GITREPO = git@github.com:kokkos/mdspan.git
include ${MAKEINCLUDES}/Make.git

