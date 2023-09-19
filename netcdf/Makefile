################################################################
####
#### Makefile for Netcdf installation
####
#### https://netcdf.geodynamics.org/
#### https://github.com/geodynamics/netcdf
####
#### IS THIS STILL VALID?
#### SUPERSEDE WITH PNETCDF?
####
################################################################

PACKAGE = NETCDF
PACKAGEVERSION = 4.9.2
MODE = mpi
MODULES = phdf5

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install :
	@make configure build JCOUNT=${JCOUNT}

TGZURL = https://downloads.unidata.ucar.edu/netcdf-c/${PACKAGEVERSION}/netcdf-c-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

