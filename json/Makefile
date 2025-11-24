################################################################
####
#### Makefile for JSON installation
####
#### Nvidia problem:
#### https://github.com/nlohmann/json/discussions/4899
####
################################################################

PACKAGE = json
ABOUT = JSON for Modern C++
PACKAGEVERSION = 3.12.0
URL = https://github.com/nlohmann/json
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

DOWNLOADURL = https://github.com/nlohmann/json/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

