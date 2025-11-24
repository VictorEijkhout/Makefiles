################################################################
####
#### Makefile for Gsl installation
####
################################################################

PACKAGE = GSL
ABOUT = GNU Scientific Library
URL = http://www.gnu.org/software/gsl
PACKAGEVERSION = 2.8
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
CONFIGUREFLAGS = \
    --with-pic \
    --enable-shared \
    --disable-dependency-tracking 

PKGCONFIGLIB = pkgconfig
HASBIN = 1
include ${MAKEINCLUDES}/Make.install
.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build

DOWNLOADURL = https://ftp.gnu.org/gnu/gsl/gsl-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

