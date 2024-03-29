################################################################
####
#### Makefile for Netcdf-fortran installation
####
################################################################

PACKAGE = NETCDF-FORTRAN
URL = https://github.com/Unidata/netcdf-fortran \
    https://docs.unidata.ucar.edu/netcdf-c/current/building_netcdf_fortran.html
PACKAGEVERSION = 4.6.1
MODE = mpi
MODULES = phdf5 pnetcdf

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
CMAKEPREFIXPATHSET = 1
DEPENDSON = netcdf
PKGCONFIGLIB = pkgconfig
NETCDF_PARALLEL = ON
CMAKEFLAGS = \
    -D netCDF_INCLUDE_DIR=${TACC_NETCDF_INC} \
    -D netCDF_LIBRARIES=${TACC_NETCDF_LIB}/libnetcdf.so \
    -D ENABLE_DAP=OFF
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: seq par
info ::
	@echo "make seq par : separate versions"
seq :
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    MODE=seq MODULES="hdf5 netcdf" NETCDF_PARALLEL=OFF \
	    MODULENAME=netcdff MODULENAMEALT=netcdf DEPENDSON=netcdf
par :
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    MODE=mpi MODULES="phdf5 pnetcdf" NETCDF_PARALLEL=ON \
	    MODULENAME=pnetcdff MODULENAMEALT=netcdf DEPENDSON=pnetcdf
.PHONY: default_install
info ::
	@echo "make default_install: parallel version"
default_install : par 

# TGZURL = https://downloads.unidata.ucar.edu/netcdf-c/${PACKAGEVERSION}/netcdf-c-${PACKAGEVERSION}.tar.gz
TGZURL = https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

