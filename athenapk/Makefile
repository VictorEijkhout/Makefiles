################################################################
####
#### Makefile for Athenapk installation
####
################################################################

PACKAGE = athenapk
URL = https://github.com/parthenon-hpc-lab/athenapk
PACKAGEVERSION = git

MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

##
## cmake version
##
CMAKEPREFIXPATHSET = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

GITREPO = https://github.com/parthenon-hpc-lab/athenapk.git
SUBMODULE = 1
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

