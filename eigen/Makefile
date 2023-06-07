################################################################
####
#### Makefile for Eigen installation
####
################################################################

PACKAGE = EIGEN
PACKAGEVERSION = 3.4.0
MODULES = 
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIG = share/pkgconfig
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://gitlab.com/libeigen/eigen/-/archive/${PACKAGEVERSION}/eigen-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download



