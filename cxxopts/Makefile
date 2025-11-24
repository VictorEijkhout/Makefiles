################################################################
####
#### Makefile for CXXOPTS installation
####
################################################################

PACKAGE = CXXOPTS
ABOUT = Lightweight C++ command line option parser
PACKAGEVERSION = 3.3.1
URL = https://github.com/jarro2783/cxxopts
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
NOLIB = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build public

DOWNLOADURL = https://github.com/jarro2783/cxxopts/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

####
#### gcc 15 problem
#### https://github.com/jarro2783/cxxopts/issues/455
####
BOOTSTRAP = sed -i -e 32i\"\#include <cstdint>\"  include/cxxopts.hpp

GITREPO = https://github.com/jarro2783/cxxopts.git
include ${MAKEINCLUDES}/Make.git
