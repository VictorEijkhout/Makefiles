################################################################
####
#### Makefile for SZ installation
####
################################################################

PACKAGE = SZ
PACKAGEVERSION = 2.1.12
PACKAGEVERSIONDOWNLOAD = 2.1.12.5
URL = https://github.com/szcompressor/SZ/
MODE = seq

include ${MAKEINCLUDES}/Make.info

CMAKFLAGS=
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
PKGCONFIG = share/pkgconfig
CMAKEPREFIXPATHSET = 1
POSTINSTALLinLIB = ln libSZ.so libsz.so

default_install : configure build

##
## Download and clone
##
TGZURL = https://github.com/szcompressor/SZ/archive/refs/tags/v${PACKAGEVERSIONDOWNLOAD}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/szcompressor/SZ.git
include ${MAKEINCLUDES}/Make.git

##
## Clean
##
include ${MAKEINCLUDES}/Make.clean
