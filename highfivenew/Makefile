################################################################
####
#### Makefile for Highfive installation
####
################################################################

PACKAGE = highfive
ABOUT = Header-only C++ HDF5 interface
PACKAGEVERSION = 3.2.0
# DISCONTINUED 2.10.1
URL = https://github.com/highfive-devs/highfive
# DISCONTINUED https://github.com/BlueBrain/HighFive
DOCURL = https://highfive-devs.github.io/highfive/
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
    -D HIGHFIVE_UNIT_TESTS=OFF \
    -D HIGHFIVE_PARALLEL_HDF5=${HDF5PARALLEL} \
    -D HIGHFIVE_USE_BOOST=ON \
    -D HIGHFIVE_USE_EIGEN=ON \
    -D HIGHFIVE_USE_XTENSOR=OFF \
    -D HIGHFIVE_USE_HALF_FLOAT=OFF
ifeq "${TACC_FAMILY_COMPILER}" "gcc"
  CMAKEFLAGS += -DCMAKE_EXE_LINKER_FLAGS='-lpthread -lgfortran'
else
  CMAKEFLAGS += -DCMAKE_EXE_LINKER_FLAGS='-lpthread'
endif
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
	    HDF5PARALLEL=ON  MODE=mpi MODULES="${MODULES_par}" \
	    MODULENAME=phighfive MODULENAMEALT=highfive
info ::
	@echo "make default_install : sequential version"
default_install :: seq
#
# module stuff
#
DEPENDSONCURRENT=hdf5
DEPENDSON = boost eigen

#
# Download
#
DOWNLOADURL = https://github.com/highfive-devs/highfive/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
# DISCONTINUED https://github.com/BlueBrain/HighFive/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/highfive-devs/highfive.git
SUBMODULE = 1
include ${MAKEINCLUDES}/Make.git
