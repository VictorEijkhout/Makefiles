################################################################
####
#### Makefile for Netcdf installation
####
#### this is the C version
####
################################################################

PACKAGE = NETCDF
URL = https://netcdf.geodynamics.org/ \
    https://github.com/Unidata/netcdf \
    https://docs.unidata.ucar.edu/nug/current/getting_and_building_netcdf.html#configure_options
PACKAGEVERSION = 4.9.2
MODE = mpi
MODULES = phdf5

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
CMAKEPREFIXPATHSET = 1
PKGCONFIGLIB = pkgconfig
NETCDF_PARALLEL = ON
DEPENDSON = phdf5

##
## NOTE the enable_pnetcdf flag is for parallel-netcdf 1.9
##
CMAKEFLAGS = \
    -D ENABLE_PNETCDF=OFF \
    -D ENABLE_PARALLEL4=${NETCDF_PARALLEL} \
    -D ENABLE_PARALLEL_TESTS=OFF \
    -D ENABLE_NETCDF_4=ON \
    -D ENABLE_DAP=OFF
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: seq par
info ::
	@echo "make seq par : separate versions"
seq :
	@if [ ! -z "${TACC_PHDF5_DIR}" ] ; then \
	    echo "Load sequential hdf5 instead of parallel" && exit 1 ; fi
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    MODE=seq MODULES=hdf5 NETCDF_PARALLEL=OFF \
	    MODULENAME=netcdf DEPENDSON=hdf5
par :
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    MODE=mpi MODULES=phdf5 NETCDF_PARALLEL=ON \
	    MODULENAME=pnetcdf MODULENAMEALT=netcdf DEPENDSON=phdf5
.PHONY: default_install
info ::
	@echo "make default_install: parallel version, also available: seq"
default_install : par

TGZURL = https://downloads.unidata.ucar.edu/netcdf-c/${PACKAGEVERSION}/netcdf-c-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = http://github.com/Unidata/netcdf-c
include ${MAKEINCLUDES}/Make.git
include ${MAKEINCLUDES}/Make.clean
