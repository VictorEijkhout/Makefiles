################################################################
####
#### Makefile for Netcdf installation
####
#### https://netcdf.geodynamics.org/
#### https://github.com/geodynamics/netcdf
####
################################################################

PACKAGE = NETCDF
PACKAGEVERSION = 4.9.2
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

notCMAKEFLAGS = 

include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://downloads.unidata.ucar.edu/netcdf-c/${PACKAGEVERSION}/netcdf-c-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

