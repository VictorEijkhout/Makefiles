################################################################
####
#### Makefile for GDAL installation
####
################################################################

PACKAGE = GDAL
URL = https://gdal.org/index.html
PACKAGEVERSION = 3.7.0
MODE = seq
MODULES = jsonc proj swig zlib

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
HASBIN = 1
# the following solves:
# libsqlite3.so lacks the RTree extension!
CMAKEFLAGS = \
    -D CMAKE_BUILD_TYPE=Debug \
    \
    -D GDAL_BUILD_OPTIONAL_DRIVERS=OFF \
    -D OGR_BUILD_OPTIONAL_DRIVERS=OFF \
    -D GDAL_USE_CURL=OFF \
    -D GDAL_USE_CRYPTOPP=OFF \
    -D GDAL_USE_OPENSSL=OFF \
    \
    -D ACCEPT_MISSING_SQLITE3_RTREE:BOOL=ON \
    -D GDAL_USE_JSONC_INTERNAL=OFF
include ${MAKEINCLUDES}/Make.cmake
EXTRAINSTALLVARS = GDAL_DATA=share/gdal
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install : configure build public

TGZURL = https://github.com/OSGeo/gdal/releases/download/v${PACKAGEVERSION}/gdal-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
untar ::
	@source ${MAKEINCLUDES}/names.sh \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEB\ASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 && f=ogrgeopackagetablelayer.cpp \
	 && echo " .. fixing source file $$f" \
	 && cd $${srcdir}/ogr/ogrsf_frmts/gpkg/ \
	 && mv $$f tmp \
	 && cat $${scriptdir}/gdal_fix.txt tmp > $$f
