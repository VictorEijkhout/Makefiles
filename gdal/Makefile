################################################################
####
#### Makefile for GDAL installation
#### https://gdal.org/index.html
####
################################################################

PACKAGE = GDAL
PACKAGEVERSION = 3.7.0
MODE = seq
MODULES = jsonc proj zlib

include ${MAKEINCLUDES}/Make.info

# ifeq "${TACC_FAMILY_COMPILER}" "intel"
#   LIBDIR = lib64
#   PKGCONFIG = lib64/pkgconfig
# else
#   LIBDIR = lib
#   PKGCONFIG = lib/pkgconfig
# endif
HASBIN = 1
# the following solves:
# libsqlite3.so lacks the RTree extension!
CMAKEFLAGS = \
    -D CMAKE_BUILD_TYPE=Release \
    -D ACCEPT_MISSING_SQLITE3_RTREE:BOOL=ON \
    -D GDAL_USE_JSONC_INTERNAL=OFF
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/OSGeo/gdal/releases/download/v3.7.0/gdal-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
