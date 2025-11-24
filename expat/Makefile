################################################################
####
#### Makefile for EXPAT installation
####
################################################################

PACKAGE = EXPAT
ABOUT = stream-oriented XML parser library written in C.
PACKAGEVERSION = 2.6.2
PACKAGEVERSIONDOWNLOAD = 2_6_2
URL = https://libexpat.github.io/
MODE = seq

include ${MAKEINCLUDES}/Make.info

##
## Installation
##
## for clones only: BEFORECONFIGURECMDS = ./buildconf.sh
include ${MAKEINCLUDES}/Make.configure
PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
HASBIN = 1
include ${MAKEINCLUDES}/Make.install
.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build public

##
## Download
##
DOWNLOADURL = https://github.com/libexpat/libexpat/releases/download/R_${PACKAGEVERSIONDOWNLOAD}/expat-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

##
## Clean
##
include ${MAKEINCLUDES}/Make.clean

