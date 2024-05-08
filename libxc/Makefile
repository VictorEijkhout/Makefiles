################################################################
####
#### Makefile for libxc
####
#### Note: Conquest only supports version 4 and 5, but version 6 is out
####
################################################################

PACKAGE = LIBXC
PACKAGEVERSION = 6.2.2
URL = https://gitlab.com/libxc/libxc
# 5.2.3
MODE = seq

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
    -D CMAKE_C_FLAGS="-O2" \
    -D ENABLE_FORTRAN=ON \
    -D ENABLE_XHOST=OFF \
    -D BUILD_SHARED_LIBS=OFF
## LIBDIR = lib64
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1

TGZURL = https://gitlab.com/libxc/libxc/-/archive/${PACKAGEVERSION}/libxc-${PACKAGEVERSION}.tar.gz
## http://www.tddft.org/programs/libxc/down.php?file=${PACKAGEVERSION}/libxc-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
