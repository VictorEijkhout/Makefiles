################################################################
####
#### Makefile for EXPAT installation
####
################################################################

PACKAGE = EXPAT
PACKAGEVERSION = 2.6.2
PACKAGEVERSIONDOWNLOAD = 2_6_2
URL = https://libexpat.github.io/
MODE = seq

include ${MAKEINCLUDES}/Make.info

##
## Installation
##
## BEFORECONFIGURECMDS = ./buildconf.sh
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build public

##
## Download
##
TGZURL = https://github.com/libexpat/libexpat/releases/download/R_${PACKAGEVERSIONDOWNLOAD}/expat-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

##
## Clean
##
include ${MAKEINCLUDES}/Make.clean

