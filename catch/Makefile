################################################################
####
#### Makefile for Catch2 installation
####
################################################################

PACKAGE = CATCH2
URL = https://github.com/catchorg/Catch2
PACKAGEVERSION = 3.5.4
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

PKGCONFIG = share/pkgconfig
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
default_install : configure build 

TGZURL = https://github.com/catchorg/Catch2/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/catchorg/Catch2.git
include ${MAKEINCLUDES}/Make.git
