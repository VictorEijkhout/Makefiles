################################################################
####
#### Makefile for Netcdf installation
####
#### IS THIS STILL VALID?
#### SUPERSEDE WITH PNETCDF?
####
################################################################

PACKAGE = NETCDF
URL = https://netcdf.geodynamics.org/ \
    https://github.com/geodynamics/netcdf \
    https://docs.unidata.ucar.edu/nug/current/getting_and_building_netcdf.html#configure_options
PACKAGEVERSION = 4.9.2
MODE = mpi
MODULES = phdf5

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
seq :
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    MODE=seq MODULES=hdf5 NETCDF_PARALLEL=OFF \
	    MODULENAME=netcdf
par :
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    MODE=mpi MODULES=phdf5 NETCDF_PARALLEL=ON \
	    MODULENAME=pnetcdf MODULENAMEALT=netcdf
.PHONY: default_install
default_install : par

TGZURL = https://downloads.unidata.ucar.edu/netcdf-c/${PACKAGEVERSION}/netcdf-c-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

