################################################################
####
#### Makefile for Mapl installation
####
#### https://stackoverflow.com/questions/45177196/how-to-use-cmake-module-path-correctly-on-windows
####
################################################################

PACKAGE = mapl
ABOUT = Foundation layer of the GEOS architecture
URL = https://github.com/GEOS-ESM/MAPL
PACKAGEVERSION = 2.56.0
# 2.64.0 is out

MODE = mpi
MODULES = gftl gftlshared pflogger yafyaml \
    ecbuild

include ${MAKEINCLUDES}/Make.info

##
## cmake version
##
include ${MAKEINCLUDES}/Make.cmake
CMAKENAME = ecbuild
include ${MAKEINCLUDES}/Make.cbuild

.PHONY: default_install
default_install : configure build

DOWNLOADURL = https://github.com/GEOS-ESM/MAPL/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean

