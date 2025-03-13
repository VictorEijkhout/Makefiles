################################################################
####
#### Makefile for libxc
####
#### Note: Conquest only supports version 4 and 5, but version 6 is out
####
################################################################

PACKAGE = LIBXC
PACKAGEVERSION = 7.0.0
# 6.2.2
# 5.2.3
URL = https://gitlab.com/libxc/libxc
MODE = seq

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
    -D CMAKE_C_FLAGS="-O2" \
    -D ENABLE_FORTRAN=ON \
    -D ENABLE_XHOST=OFF 
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
#BUILDSTATICLIBS=1
HASBIN = 1
PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

TGZURL = https://gitlab.com/libxc/libxc/-/archive/${PACKAGEVERSION}/libxc-${PACKAGEVERSION}.tar.gz
## http://www.tddft.org/programs/libxc/down.php?file=${PACKAGEVERSION}/libxc-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
