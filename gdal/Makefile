################################################################
####
#### Makefile for GDAL installation
####
################################################################

PACKAGE = GDAL
PACKAGEVERSION = 3.7.0
MODE = seq

include ${MAKEINCLUDES}/Make.info

# ifeq "${TACC_FAMILY_COMPILER}" "intel"
#   LIBDIR = lib64
#   PKGCONFIG = lib64/pkgconfig
# else
#   LIBDIR = lib
#   PKGCONFIG = lib/pkgconfig
# endif
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/OSGeo/gdal/releases/download/v3.7.0/gdal-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
