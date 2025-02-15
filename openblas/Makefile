################################################################
####
#### Makefile for OPENBLAS installation
####
################################################################

PACKAGE = OPENBLAS
URL = https://github.com/OpenMathLib/OpenBLAS
PACKAGEVERSION = 0.3.29
MODE = seq

include ${MAKEINCLUDES}/Make.info

## PKGCONFIGLIB = pkgconfig

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install : configure build

TGZURL = https://github.com/OpenMathLib/OpenBLAS/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
