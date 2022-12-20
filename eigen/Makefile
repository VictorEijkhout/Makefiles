################################################################
####
#### Makefile for Eigen installation
####
################################################################

PACKAGE = EIGEN
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 3.4.0
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODULEROOT = ${HOME}/modulefiles
MODULES = 
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIGSET = share/pkgconfig
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://gitlab.com/libeigen/eigen/-/archive/${PACKAGEVERSION}/eigen-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download



