################################################################
####
#### Makefile for Netcdf installation
####
#### this is the C version
####
#### pnetcdf handling: https://github.com/Unidata/netcdf-c/issues/2872
####
################################################################

PACKAGE = NETCDF
URL = https://netcdf.geodynamics.org/ \
    https://github.com/Unidata/netcdf-c \
    https://docs.unidata.ucar.edu/nug/current/getting_and_building_netcdf.html#configure_options
PACKAGEVERSION = 4.9.2
# 4.9.3 is out
MODE = seq
MODULES_seq = hdf5
MODULES_par = pnetcdf phdf5
MODULES = ${MODULES_seq}

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
CMAKEPREFIXPATHSET = 1
PKGCONFIGLIB = pkgconfig
NETCDF_PARALLEL = OFF
DEPENDSON = hdf5
LINKLIB64toLIB = 1
CPTOINSTALLDIR = examples
##
## NOTE the enable_pnetcdf flag is for parallel-netcdf 1.9
##
CMAKEFLAGS = \
    -D NOT_ENABLE_PNETCDF=${NETCDF_PARALLEL} ${CMAKE_C_EXTRA_FLAGS} \
    -D ENABLE_PARALLEL4=${NETCDF_PARALLEL} \
    -D ENABLE_PARALLEL_TESTS=OFF \
    -D ENABLE_NETCDF_4=ON \
    -D ENABLE_DAP=ON
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: seq par
info ::
	@echo "make seq par : separate versions"
seq :
	@if [ ! -z "${TACC_PHDF5_DIR}" ] ; then \
	    echo "Load sequential hdf5 instead of parallel" && exit 1 ; fi
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    MODE=seq MODULES="${MODULES_seq}" NETCDF_PARALLEL=OFF \
	    MODULENAME=netcdf DEPENDSON=hdf5
par :
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    MODE=mpi MODULES="${MODULES_par}" NETCDF_PARALLEL=ON \
	    CMAKE_C_EXTRA_FLAGS="-D CMAKE_C_FLAGS='-I${TACC_PHDF5_INC} -I${TACC_PNETCDF_INC}' -D CMAKE_LINKER_FLAGS='-L${TACC_PNETCDF_LIB} -lpnetcdf'" \
	    MODULENAME=parallel-netcdf MODULEVARNAME=parallelnetcdf DEPENDSON=phdf5
.PHONY: default_install
info ::
	@echo "make default_install: sequential version, also available: par"
default_install : seq

TGZURL = https://downloads.unidata.ucar.edu/netcdf-c/${PACKAGEVERSION}/netcdf-c-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = http://github.com/Unidata/netcdf-c
include ${MAKEINCLUDES}/Make.git
include ${MAKEINCLUDES}/Make.clean
