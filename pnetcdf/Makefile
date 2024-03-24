################################################################
####
#### Makefile for Netcdf installation
####
#### https://parallel-netcdf.github.io/
####
################################################################

PACKAGE = pnetcdf
URL = https://parallel-netcdf.github.io/ \
    https://github.com/Parallel-NetCDF/PnetCDF
PACKAGEVERSION = 1.12.3
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

GITREPO = https://github.com/Parallel-NetCDF/PnetCDF.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

