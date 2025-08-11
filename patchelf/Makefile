################################################################
####
#### Makefile for PATCHELF installation
####
################################################################

PACKAGE = PATCHELF
ABOUT = A small utility to modify the dynamic linker and RPATH of ELF executables
PACKAGEVERSION = 0.18.0
URL = https://github.com/NixOS/patchelf
MODE = core

include ${MAKEINCLUDES}/Make.info

NOLIB = 1
NOINCLUDE = 1
HASBIN = 1
BEFORECONFIGURECMDS = ./bootstrap.sh
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build public

TGZURL = https://github.com/NixOS/patchelf/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/jarro2783/patchelf.git
