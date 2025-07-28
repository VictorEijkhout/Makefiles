################################################################
####
#### Makefile for Openmm installation
####
#### internal compiler error:
#### https://github.com/openmm/openmm/issues/4843
####
################################################################

PACKAGE = openmmcuda
URL = https://github.com/openmm/openmm
PACKAGEVERSION = 8.2.0
# 8.3.1 is out

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
    -D OPENMM_BUILD_PYTHON_WRAPPERS=ON \
    -D OPENMM_BUILD_C_AND_FORTRAN_WRAPPERS=ON
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

TGZURL = https://github.com/openmm/openmm/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = git@github.com:openmm/openmm.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean
