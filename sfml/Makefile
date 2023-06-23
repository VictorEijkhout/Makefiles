################################################################
####
#### Makefile for SFML installation
####
################################################################

PACKAGE = SFML
URL = https://github.com/SFML/SFML
PACKAGEVERSION = 2.6.0
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/SFML/SFML/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/SFML/SFML.git
include ${MAKEINCLUDES}/Make.git
