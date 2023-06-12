################################################################
####
#### Makefile for Netcdf installation
####
#### https://parallel-netcdf.github.io/
#### https://github.com/Parallel-NetCDF/PnetCDF
####
################################################################

PACKAGE = pNETCDF
PACKAGEVERSION = 1.12.3
MODE = mpi
MODULES = hdf5

include ${MAKEINCLUDES}/Make.info
CONFIGUREFLAGS = --enable-fortran --enable-shared
include ${MAKEINCLUDES}/Make.configure
HASBIN=1
include ${MAKEINCLUDES}/Make.install

TGZURL = https://parallel-netcdf.github.io/Release/pnetcdf-${PACKAGEVERSION}.tar.gz
## this one is missing some script:
## https://github.com/Parallel-NetCDF/PnetCDF/archive/refs/tags/checkpoint.${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

