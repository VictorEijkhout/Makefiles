################################################################
####
#### Makefile for Highfive installation
####
################################################################

PACKAGE = highfive
PACKAGEVERSION = 2.7.1
URL = https://github.com/BlueBrain/HighFive
MODE = mpi
## note: modules needs hdf5 or phdf5 depending
MODULES = boost eigen 

include ${MAKEINCLUDES}/Make.info

##
## cmake version
##
HDF5PARALLEL = ON
CMAKEFLAGS = \
    -D HIGHFIVE_PARALLEL_HDF5=${HDF5PARALLEL} \
    -D HIGHFIVE_USE_BOOST=ON \
    -D HIGHFIVE_USE_EIGEN=ON \
    -D HIGHFIVE_USE_XTENSOR=OFF \
    -D HIGHFIVE_USE_HALF_FLOAT=OFF
include ${MAKEINCLUDES}/Make.cmake
## .cmake files in share/HighFive/CMake/
NOLIB = 1
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make par : parallel installation"
	@echo "make seq : sequential installation"
.PHONY: seq par
seq :
	make configure build JCOUNT=${JCOUNT} \
	    PACKAGEVERSION=${PACKAGEVERSION} \
	    HDF5PARALLEL=OFF MODE=seq MODULES="${MODULES} hdf5" MODULENAME=highfive
par :
	make configure build JCOUNT=${JCOUNT} \
	    PACKAGEVERSION=${PACKAGEVERSION} \
	    HDF5PARALLEL=ON  MODE=mpi MODULES="${MODULES} phdf5" MODULENAME=phighfive

TGZURL = https://github.com/BlueBrain/HighFive/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = git@github.com:BlueBrain/HighFive.git
include ${MAKEINCLUDES}/Make.git
