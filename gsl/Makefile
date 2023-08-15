################################################################
####
#### Makefile for GSL-LITE installation
####
################################################################

PACKAGE = MSgsl
PACKAGEVERSION = 4.0.0
URL = https://github.com/microsoft/gsl
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
NOLIB = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/microsoft/GSL/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
