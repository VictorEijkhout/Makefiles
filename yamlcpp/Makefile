################################################################
####
#### Makefile for YAMLCPP installation
####
################################################################

PACKAGE = YAMLCPP
ABOUT = A YAML parser and emitter in C++ 
PACKAGEVERSION = 0.9.0
URL = https://github.com/jbeder/yaml-cpp
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build public

####
#### Download
####
DOWNLOADURL = https://github.com/jbeder/yaml-cpp/archive/refs/tags/yaml-cpp-${PACKAGEVERSION}.tar.gz
# https://github.com/jbeder/yaml-cpp/issues/1359
# should be fixed in 0.9.0
# BOOTSTRAP = sed -i -e 3i\"\#include <cstdint>\"  src/*.cpp
include ${MAKEINCLUDES}/Make.download
