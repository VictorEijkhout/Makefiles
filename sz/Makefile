################################################################
####
#### Makefile for SZ installation
####
################################################################

PACKAGE = SZ
PACKAGEVERSION = 2.1.12
PACKAGEVERSIONDOWNLOAD = 2.1.12.5
MODE = seq

include ${MAKEINCLUDES}/Make.info

CMAKFLAGS=
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
POSTINSTALLinLIB = ln libSZ.so libsz.so

default_install : configure build

TGZURL = https://github.com/szcompressor/SZ/archive/refs/tags/v${PACKAGEVERSIONDOWNLOAD}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

GITREPO = https://github.com/szcompressor/SZ.git
