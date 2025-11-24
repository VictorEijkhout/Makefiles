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
PKGCONFIGLIB = pkgconfig
LINKLIB64toLIB = 1
NETCDF_PARALLEL = OFF
CPTOINSTALLDIR = examples
CMAKEFLAGS_seq = \
    -D netCDF_INCLUDE_DIR=${TACC_NETCDF_INC} \
    -D netCDF_LIBRARIES=${TACC_NETCDF_LIB}/libnetcdf.so \
    -D ENABLE_DAP=OFF
CMAKEFLAGS_par = \
    -D netCDF_INCLUDE_DIR=${TACC_PARALLELNETCDF_INC} \
    -D netCDF_LIBRARIES=${TACC_PARALLELNETCDF_LIB}/libnetcdf.so \
    -D ENABLE_DAP=OFF
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: seq par
info ::
	@echo "make seq par : separate versions"
seq :
	@make --no-print-directory configure build \
	    NETCDF_PARALLEL=OFF JCOUNT=${JCOUNT} \
	    CMAKEFLAGS="${CMAKEFLAGS_$@}" \
	    MODE=seq MODULES="hdf5 netcdf" \
	    DEPENDSON=netcdf \
	    MODULENAME=netcdff  
par :
	@make --no-print-directory configure build \
	    NETCDF_PARALLEL=ON JCOUNT=${JCOUNT} \
	    CMAKEFLAGS="${CMAKEFLAGS_$@}" \
	    MODE=mpi MODULES="parallel-netcdf phdf5" \
	    MODULENAME=parallel-netcdff MODULEVARNAME=PARALLELNETCDFF

.PHONY: default_install
info ::
	@echo "make default_install: sequential version"
default_install : seq

# DOWNLOADURL = https://downloads.unidata.ucar.edu/netcdf-c/${PACKAGEVERSION}/netcdf-c-${PACKAGEVERSION}.tar.gz
DOWNLOADURL = https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

