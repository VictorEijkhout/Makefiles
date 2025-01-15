################################################################
####
#### Makefile for Blis installation
####
################################################################

PACKAGE = blis
PACKAGEVERSION = 1.0
# git
URL = https://github.com/flame/blis
MODE = seq

include ${MAKEINCLUDES}/Make.info

CONFIGUREFLAGS = generic
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
PKGCONFIG = share/pkgconfig
CPTOINSTALLDIR = examples

.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build public

TGZURL = https://github.com/flame/blis/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/flame/blis.git
include ${MAKEINCLUDES}/Make.git

