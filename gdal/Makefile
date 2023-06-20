################################################################
####
#### Makefile for GDAL installation
####
################################################################

PACKAGE = GDAL
URL = https://gdal.org/index.html
PACKAGEVERSION = 3.7.0
MODE = seq
MODULES = jsonc proj zlib

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
HASBIN = 1
# the following solves:
# libsqlite3.so lacks the RTree extension!
CMAKEFLAGS = \
    -D CMAKE_BUILD_TYPE=Release \
    -D ACCEPT_MISSING_SQLITE3_RTREE:BOOL=ON \
    -D GDAL_USE_JSONC_INTERNAL=OFF
include ${MAKEINCLUDES}/Make.cmake
EXTRAINSTALLVARS = GDAL_DATA=share/gdal
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/OSGeo/gdal/releases/download/v${PACKAGEVERSION}/gdal-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
