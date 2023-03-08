################################################################
####
#### Makefile for Libisl installation
####
#### https://libisl.sourceforge.io/
####
################################################################

PACKAGE = libisl
PACKAGEVERSION = 0.25
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIGSET = lib/pkgconfig
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install

TGZURL = https://libisl.sourceforge.io/isl-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

