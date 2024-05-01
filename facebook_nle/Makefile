################################################################
####
#### Makefile for NLE installation
####
################################################################

PACKAGE = NLE
PACKAGEVERSION = 0.9.0
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
# NOLIB = 1

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
default_install : configure build public

TGZURL = https://github.com/facebookresearch/nle/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
