################################################################
####
#### Makefile for PROJ installation
#### https://proj.org/index.html
####
################################################################

PACKAGE = PROJ
URL = https://proj.org/en/9.3/download.html
PACKAGEVERSION = 9.3.1
MODE = seq
MODULES = sqlite

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = --trace \
    -D BUILD_TESTING=OFF \
    -D EXE_SQLITE3=${TACC_SQLITE_BIN}/sqlite3 \
    -D SQLITE3_INCLUDE_DIR=${TACC_SQLITE_INC} \
    -D SQLITE3_LIBRARY=${TACC_SQLITE_LIB}/libsqlite3.so

CMAKEPREFIXPATHSET = 1
PKGCONFIGLIB = pkgconfig
HASBIN = 1
DEPENDSON = sqlite
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install : configure build public

TGZURL = https://download.osgeo.org/proj/proj-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/OSGeo/PROJ.git
include ${MAKEINCLUDES}/Make.git
