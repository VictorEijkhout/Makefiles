################################################################
####
#### Makefile for Netcdf installation
####
################################################################

PACKAGE = pnetcdf
URL = https://parallel-netcdf.github.io/ \
    https://github.com/Parallel-NetCDF/PnetCDF
PACKAGEVERSION = 1.14.0
MODE = mpi
MODULES = automake

include ${MAKEINCLUDES}/Make.info
CONFIGUREFLAGS = --enable-fortran --enable-shared
FORCERECONF = 1
# WEIRD only for gcc15 in jail.
ifeq "${TACC_FAMILY_COMPILER}" "gcc"
  BEFORECONFIGURECMDS = export CFLAGS=\"-std=c99 -Wno-error -Wno-incompatible-pointer-types -Wno-implicit-function-declaration -Wno-int-conversion\"
endif
include ${MAKEINCLUDES}/Make.configure
HASBIN=1
PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.install

info ::
	@echo "make default_install : configure install"
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

