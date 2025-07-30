################################################################
####
#### Makefile for CMAKE installation
####
################################################################

PACKAGE = CMAKE
PACKAGEVERSION = 3.31.5
URL = https://cmake.org/
MODE = core

include ${MAKEINCLUDES}/Make.info

CMAKEPREFIXPATHSET = 1
NOLIB = 1
NOINC = 1
HASBIN = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
ifeq "${TACC_SYSTEM}" "frontera"
  LDLIBRARYPATHADDABS = /home1/apps/gcc/13.2.0/lib64
endif
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build 

TGZURL = https://github.com/Kitware/CMake/releases/download/v${PACKAGEVERSION}/cmake-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
