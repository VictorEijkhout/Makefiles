################################################################
####
#### Makefile for Zoltan installation
####
#### UNFINISHED
####
################################################################

PACKAGE = ZOLTAN
URL = https://sandialabs.github.io/Zoltan/
PACKAGEVERSION = git
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
# HASBIN=1
PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure install"
default_install ::  configure install

##
## download
##
GITREPO = https://github.com/sandialabs/Zoltan.git
include ${MAKEINCLUDES}/Make.git

