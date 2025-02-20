################################################################
####
#### Makefile for SuiteSparse installation
####
################################################################

PACKAGE = SUITESPARSE
PACKAGEVERSION = 7.8.3
URL = https://github.com/DrTimothyAldenDavis/SuiteSparse
MODE = seq
MODULES = mpfr

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

CMAKEFLAGS=\
        -D CMAKE_C_COMPILER=$${cc} -D CMAKE_CXX_COMPILER=$${cxx} \
        -D CMAKE_C_FLAGS=-diag-disable=10441 -D CMAKE_CXX_FLAGS=-diag-disable=10441
HASBIN = 
CMAKEPREFIXPATHSET = 1
.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build
DEPENDSON = gmp,mpfr

TGZURL = https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

