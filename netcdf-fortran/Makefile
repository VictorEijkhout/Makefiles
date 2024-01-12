################################################################
####
#### Makefile for Netcdf-fortran installation
####
################################################################

PACKAGE = NETCDF-FORTRAN
URL = https://netcdf.geodynamics.org/ \
    https://github.com/geodynamics/netcdf \
    https://docs.unidata.ucar.edu/nug/current/getting_and_building_netcdf.html#configure_options
PACKAGEVERSION = 4.6.1
MODE = mpi
MODULES = phdf5 pnetcdf

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
PKGCONFIGLIB = pkgconfig
NETCDF_PARALLEL = ON
CMAKEFLAGS = \
    -D ENABLE_PARALLEL4=${NETCDF_PARALLEL} \
    -D ENABLE_NETCDF_4=ON \
    -D ENABLE_DAP=OFF
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: seq par
info ::
	@echo "make seq par : separate versions"
seq :
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    MODE=seq MODULES="hdf5 netcdf" NETCDF_PARALLEL=OFF \
	    MODULENAME=netcdf-fortran
par :
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    MODE=mpi MODULES="phdf5 pnetcdf" NETCDF_PARALLEL=ON \
	    MODULENAME=pnetcdf-fortran MODULENAMEALT=netcdf
.PHONY: default_install
info ::
	@echo "make default_install: parallel version"
default_install : par

# TGZURL = https://downloads.unidata.ucar.edu/netcdf-c/${PACKAGEVERSION}/netcdf-c-${PACKAGEVERSION}.tar.gz
TGZURL = https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

