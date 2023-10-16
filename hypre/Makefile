################################################################
####
#### Makefile for Hypre installation
####
################################################################

PACKAGE = hypre
URL = https://github.com/hypre-space/hypre
PACKAGEVERSION = 2.29.0

MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "make small big ( 32/64 bit int versions)"
.PHONY: small big
small :
	@make --no-print-directory configure build public JCOUNT=${JCOUNT} \
	    INSTALLVARIANT=i32 
big :
	@make --no-print-directory configure build public JCOUNT=${JCOUNT} \
	    CONFIGUREFLAGS=--enable-bigint \
	    INSTALLVARIANT=i64 INSTALLEXT=i64
default_install : small big

##
## configure version
##
CONFIGURESUBDIR = src
MAKESUBDIR = src
#include ${MAKEINCLUDES}/Make.configure
#include ${MAKEINCLUDES}/Make.install

##
## cmake version
##
CMAKESUBDIR = src
CMAKEFLAGS = -D HYPRE_ENABLE_SHARED=ON
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/hypre-space/hypre/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean

