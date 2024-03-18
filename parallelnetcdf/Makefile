################################################################
####
#### Makefile for Parallel-Netcdf installation
####
################################################################

PACKAGE = PARALLELNETCDF
URL = https://parallel-netcdf.github.io/ \
    https://github.com/Parallel-NetCDF/PnetCDF
PACKAGEVERSION = 1.12.3
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

AUTORECONF = autoreconf --install
include ${MAKEINCLUDES}/Make.configure

HASBIN=1
CONFIGUREFLAGS="--enable-shared"
PKGCONFIGLIB = pkgconfig
CPTOINSTALLDIR = examples
include ${MAKEINCLUDES}/Make.install
.PHONY: default_install
info ::
	@echo "make default_install: parallel version, also available: seq"
default_install : configure build

TGZURL = https://parallel-netcdf.github.io/Release/pnetcdf-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/Parallel-NetCDF/PnetCDF.git
include ${MAKEINCLUDES}/Make.git
include ${MAKEINCLUDES}/Make.clean
