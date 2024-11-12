################################################################
####
#### Makefile for Netcdf-fortran installation
####
################################################################

##
## note that we have a separate name here
## but as a system rpm we will merge C & F
##
PACKAGE = netcdff
URL = https://github.com/Unidata/netcdf-fortran \
    https://docs.unidata.ucar.edu/netcdf-c/current/building_netcdf_fortran.html
PACKAGEVERSION = 4.6.1
MODE = seq
MODULES = hdf5 netcdf

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
CMAKEPREFIXPATHSET = 1
DEPENDSON = netcdf
PKGCONFIGLIB = pkgconfig
LINKLIB64toLIB = 1
NETCDF_PARALLEL = OFF
CPTOINSTALLDIR = examples
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
	    MODULENAME=netcdff  DEPENDSON=netcdf
par :
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    MODE=mpi MODULES="" NETCDF_PARALLEL=ON \
	    MODULENAME=parallelnetcdff

.PHONY: default_install
info ::
	@echo "make default_install: sequential version"
default_install : seq

# TGZURL = https://downloads.unidata.ucar.edu/netcdf-c/${PACKAGEVERSION}/netcdf-c-${PACKAGEVERSION}.tar.gz
TGZURL = https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

