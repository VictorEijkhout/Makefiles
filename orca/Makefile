################################################################
####
#### Makefile for Orca installation
####
#### WE ARE NOT ACTUALLY INSTALLING THIS
#### AND THIS MAKEFILE IS COMPLETELY WRONG.
#### I think I once registered on their website and then never used it
####
################################################################

PACKAGE = orcacuda
SKIPPACKAGE = 1
ABOUT = 
URL = 
PACKAGEVERSION = 6.1.1

## junk from here on down

MODE = seq
BASICMODULES = python swig 
MODULES = ${BASICMODULES} cuda

include ${MAKEINCLUDES}/Make.info

##
## cmake 
##
USEGPU = YES
CUDAFLAGS = -allow-unsupported-compiler -Xcompiler -std=gnu++14
CMAKEFLAGS = -Wno-dev \
    -D PYTHON_EXECUTABLE=${TACC_PYTHON_BIN}/python3 \
    -D ORCA_BUILD_PYTHON_WRAPPERS=ON \
    -D ORCA_BUILD_C_AND_FORTRAN_WRAPPERS=ON
HASBIN = 1
EXTRAINSTALLTARGET = PythonInstall

include ${MAKEINCLUDES}/Make.cmake
info ::
	@echo "make config-help"
config-help :
	@cd ${PACKAGEROOT}/${PACKAGE}/${PACKAGE}-${PACKAGEVERSION} \
	 && open docs-source//usersguide/library/02_compiling.rst
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : gpu version"
.PHONY: default_install gpu
default_install gpu :
	@export CUDAFLAGS="${CUDAFLAGS}" \
	 && export NVCC_APPEND_FLAGS="${CUDAFLAGS}" \
	 && make --no-print-directory configure build JCOUNT=${JCOUNT}

DOWNLOADURL = https://github.com/orca/orca/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = git@github.com:orca/orca.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean
