################################################################
####
#### Makefile for Libisl installation
####
################################################################

PACKAGE = libisl
PACKAGEVERSION = 0.27
URL = https://libisl.sourceforge.io/
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIG = lib/pkgconfig
include ${MAKEINCLUDES}/Make.configure
HASBIN =
include ${MAKEINCLUDES}/Make.install

.PHONY: default_install
info ::
	echo "make default_install : configure build public"
default_install : configure build public

TGZURL = https://libisl.sourceforge.io/isl-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

