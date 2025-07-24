################################################################
####
#### Makefile for Highfive installation
####
################################################################

PACKAGE = highfive
PACKAGEVERSION = 2.10.1
URL = https://github.com/BlueBrain/HighFive
DOCURL = https://bluebrain.github.io/HighFive/
MODE = mpi
## note: modules needs hdf5 or phdf5 depending
MODULES_seq = boost eigen hdf5
MODULES_par = boost eigen phdf5
MODULES = ${MODULES_seq}

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
CMAKEPREFIXPATHSET = 1
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
	    HDF5PARALLEL=OFF MODE=seq MODULENAME=highfive
par :
	make configure build JCOUNT=${JCOUNT} \
	    PACKAGEVERSION=${PACKAGEVERSION} \
	    HDF5PARALLEL=ON  MODE=mpi MODULES="${MODULES_par}" MODULENAME=phighfive
info ::
	@echo "make default_install : sequential version"
default_install :: 
	@make --no-print-directory seq par JCOUNT=${JCOUNT} PACKAGEVERSION=${PACKAGEVERSION} \
	    seq 

TGZURL = https://github.com/BlueBrain/HighFive/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = git@github.com:BlueBrain/HighFive.git
include ${MAKEINCLUDES}/Make.git
