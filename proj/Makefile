################################################################
####
#### Makefile for PROJ installation
#### https://proj.org/index.html
####
################################################################

PACKAGE = PROJ
ABOUT = Cartographic Projections and Coordinate Transformations Library
URL = https://github.com/OSGeo/PROJ
PACKAGEVERSION = 9.6.0
# 9.6.2 is out, fixes cstdint and others
# 9.7.1
MODE = seq
MODULES = sqlite

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = --trace \
    -D BUILD_TESTING=OFF \
    -D EXE_SQLITE3=${TACC_SQLITE_BIN}/sqlite3 \
    -D SQLITE3_INCLUDE_DIR=${TACC_SQLITE_INC} \
    -D SQLITE3_LIBRARY=${TACC_SQLITE_LIB}/libsqlite3.so

HASBIN = 1
DEPENDSON = sqlite
CMAKEPREFIXPATHSET = 1
PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install : configure build public

DOWNLOADURL = https://github.com/OSGeo/PROJ/archive/refs/tags/${PACKAGEVERSION}.tar.gz
### https://download.osgeo.org/proj/proj-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/OSGeo/PROJ.git
include ${MAKEINCLUDES}/Make.git
