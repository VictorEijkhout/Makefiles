################################################################
####
#### Makefile for Netcdf installation
####
################################################################

PACKAGE = pnetcdf
URL = https://parallel-netcdf.github.io/ \
    https://github.com/Parallel-NetCDF/PnetCDF
PACKAGEVERSION = 1.13.0
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info
CONFIGUREFLAGS = --enable-fortran --enable-shared
FORCERECONF = 1
HASBIN=1
PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.configure

include ${MAKEINCLUDES}/Make.install

info ::
	@echo "make default_install : configure install )"
.PHONY: default_install
default_install : configure build

TGZURL = https://parallel-netcdf.github.io/Release/pnetcdf-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

##
## trouble with git version
## https://github.com/Parallel-NetCDF/PnetCDF/issues/134#issuecomment-2657881334
## 
GITREPO = https://github.com/Parallel-NetCDF/PnetCDF.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

