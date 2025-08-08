################################################################
####
#### Makefile for Catch2 installation
####
################################################################

PACKAGE = CATCH2
ABOUT = A modern, C++-native, test framework for unit-tests, TDD and BDD
URL = https://github.com/catchorg/Catch2
PACKAGEVERSION = 3.8.1
# 3.9.1 is out
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

PKGCONFIG = share/pkgconfig
CMAKEPREFIXPATHSET = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
default_install : configure build 

TGZURL = https://github.com/catchorg/Catch2/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/catchorg/Catch2.git
include ${MAKEINCLUDES}/Make.git
