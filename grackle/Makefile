################################################################
####
#### Makefile for Grackle installation
####
################################################################

PACKAGE = GRACKLE
ABOUT = The Grackle chemistry and cooling library for astrophysical simulations and models.
URL = https://github.com/grackle-project/grackle
PACKAGEVERSION = 3.4.0
# 3.4.1 announced 
MODE = seq
MODULES = hdf5

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
CPTOINSTALLDIR = grackle_data_files
.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build

GITREPO = https://github.com/grackle-project/grackle.git
include ${MAKEINCLUDES}/Make.git

DOWNLOADURL = https://github.com/grackle-project/grackle/archive/refs/tags/grackle-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

