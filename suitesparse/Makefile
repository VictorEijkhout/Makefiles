################################################################
####
#### Makefile for SuiteSparse installation
####
################################################################

PACKAGE = SUITESPARSE
PACKAGEVERSION = 7.9.0
# 7.11.0 is out
URL = https://github.com/DrTimothyAldenDavis/SuiteSparse
MODE = seq
MODULES = gmp mpfr mkl

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

CMAKEFLAGS=
ifeq "${TACC_FAMILY_COMPILER}" "intel"
CMAKEBUILDTYPE = Debug
CMAKEFLAGS += \
    -D CMAKE_C_FLAGS='  -diag-disable=10441 -g -O2' \
    -D CMAKE_CXX_FLAGS='-diag-disable=10441 -g -O2'
endif
HASBIN = 1
CMAKEPREFIXPATHSET = 1
.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build
DEPENDSON = gmp,mpfr

TGZURL = https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
BOOTSTRAP = sed -i -e '/\-O3/s/^/\#VLE/' LAGraph/CMakeLists.txt

include ${MAKEINCLUDES}/Make.clean
