################################################################
####
#### Makefile for SFML installation
####
#### NOTE on MacOS only clang support
####
################################################################

PACKAGE = SFML
ABOUT  = Simple and Fast Multimedia Library
URL = https://github.com/SFML/SFML
PACKAGEVERSION = 3.0.1
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

##
## download or clone
##
# https://github.com/SFML/SFML/archive/refs/tags/3.0.1.tar.gz
TGZURL = https://github.com/SFML/SFML/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/SFML/SFML.git
include ${MAKEINCLUDES}/Make.git
