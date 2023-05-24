################################################################
####
#### Makefile for MPFR installation
####
################################################################

PACKAGE = MPFR
PACKAGEVERSION = 4.2.0
MODE = seq

include ${MAKEINCLUDES}/Make.info

## PKGCONFIGSET = lib64/pkgconfig
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://www.mpfr.org/mpfr-current/mpfr-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
