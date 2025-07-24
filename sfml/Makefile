################################################################
####
#### Makefile for SFML installation
####
#### NOTE on MacOS only clang support
####
################################################################

PACKAGE = SFML
URL = https://github.com/SFML/SFML
PACKAGEVERSION = 2.6.0
MODE = seq

include ${MAKEINCLUDES}/Make.info
info ::
	@echo "================" \
	 && echo "NOTE: on MacOS only clang supported" \
	 && echo "================"

PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

TGZURL = https://github.com/SFML/SFML/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/SFML/SFML.git
include ${MAKEINCLUDES}/Make.git
