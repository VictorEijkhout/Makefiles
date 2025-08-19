################################################################
####
#### Makefile for GDAL installation
####
################################################################

PACKAGE = GDAL
ABOUT = Translator library for raster and vector geospatial data formats
URL = https://gdal.org/index.html \
    https://github.com/OSGeo/GDAL
DOCURL = https://gdal.org/en/stable/development/building_from_source.html
PACKAGEVERSION = 3.11.0
MODE = seq
MODULES = jsonc netcdf proj swig zlib 

include ${MAKEINCLUDES}/Make.info

# weird: there is a lib/python & lib64/pkgconfig
PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
HASBIN = 1
# the following solves:
# libsqlite3.so lacks the RTree extension!
CMAKEFLAGS = \
    -D CMAKE_BUILD_TYPE=Debug \
    -D CMAKE_C_FLAGS=-Wno-deprecated-declarations \
    \
    -D GDAL_BUILD_OPTIONAL_DRIVERS=OFF \
    -D OGR_BUILD_OPTIONAL_DRIVERS=OFF \
    -D GDAL_USE_CURL=OFF \
    -D GDAL_USE_CRYPTOPP=OFF \
    -D GDAL_USE_ODBC=OFF \
    -D GDAL_USE_OPENSSL=OFF \
    \
    -D GDAL_USE_NETCDF=ON \
    -D GDAL_ENABLE_DRIVER_NETCDF=ON \
    -D GDAL_USE_JSONC_INTERNAL=OFF
##    -D ACCEPT_MISSING_SQLITE3_RTREE:BOOL=ON
include ${MAKEINCLUDES}/Make.cmake
EXTRAINSTALLVARS = GDAL_DATA=share/gdal
include ${MAKEINCLUDES}/Make.cbuild
DEPENDSON = jsonc,netcdf,proj,zlib

.PHONY: default_install
default_install : configure build public

TGZURL = https://github.com/OSGeo/gdal/releases/download/v${PACKAGEVERSION}/gdal-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
