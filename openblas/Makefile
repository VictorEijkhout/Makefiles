################################################################
####
#### Makefile for OPENBLAS installation
####
################################################################

PACKAGE = OPENBLAS
ABOUT = Toolkit for molecular simulation using high performance GPU code.
URL = https://github.com/OpenMathLib/OpenBLAS
PACKAGEVERSION = 0.3.31
MODE = seq

include ${MAKEINCLUDES}/Make.info

## PKGCONFIGLIB = pkgconfig

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install : configure build

DOWNLOADURL = https://github.com/OpenMathLib/OpenBLAS/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
