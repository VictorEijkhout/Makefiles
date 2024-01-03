################################################################
####
#### Makefile for Precice installation
####
#### only needed for the CHyPS CSA project
####
################################################################

PACKAGE = precice
URL = https://github.com/precice/precice
PACKAGEVERSION = 2.5.0

MODE = seq
MODULES = boost eigen

include ${MAKEINCLUDES}/Make.info

##
## cmake 
##
CMAKEFLAGS = -DBUILD_SHARED_LIBS=ON \
    -D PRECICE_PythonActions=OFF \
    -D CMAKE_BUILD_TYPE=Release \
    -D PRECICE_PETScMapping=OFF

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/precice/precice/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = git@github.com:precice/precice.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

