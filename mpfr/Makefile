################################################################
####
#### Makefile for MPFR installation
####
################################################################

PACKAGE = MPFR
PACKAGEVERSION = 4.2.0
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIGSET = lib/pkgconfig
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install

TGZURL = https://www.mpfr.org/mpfr-current/mpfr-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
