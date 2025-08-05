################################################################
####
#### Makefile for MPFR installation
####
################################################################

PACKAGE = MPFR
ABOUT = Efficient arbitrary-precision floating-point library
PACKAGEVERSION = 4.2.2
URL = https://www.mpfr.org/ \
    https://gitlab.inria.fr/mpfr/mpfr
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIG = lib/pkgconfig
CMAKEPREFIXPATHSET = 1
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

TGZURL = https://www.mpfr.org/mpfr-current/mpfr-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
