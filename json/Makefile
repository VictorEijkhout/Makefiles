################################################################
####
#### Makefile for JSON installation
####
#### https://github.com/nlohmann/json
####
################################################################

PACKAGE = json
PACKAGEVERSION = 3.11.2
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIG = share/pkgconfig
include ${MAKEINCLUDES}/Make.cmake
NOLIB = 1
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

TGZURL = https://github.com/nlohmann/json/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

