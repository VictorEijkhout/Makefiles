################################################################
####
#### Makefile for DELFT3D installation
####
################################################################

PACKAGE = DELFT3D
## URL = 
PACKAGEVERSION = 2024.03
MODE = mpi
MODULES = gdal hdf5 jsonc netcdf proj zlib 

include ${MAKEINCLUDES}/Make.info

#PKGCONFIGLIB = pkgconfig
#CMAKEPREFIXPATHSET = 1
HASBIN = 1

CMAKEFLAGS = 
CMAKESUBDIR = cmake
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install : configure build public

# TGZURL = 
# include ${MAKEINCLUDES}/Make.download
