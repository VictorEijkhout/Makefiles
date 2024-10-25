################################################################
####
#### Makefile for Vtkhdf installation
####
#### for ticket
#### https://consult.tacc.utexas.edu/Ticket/Display.html?id=102066
####
################################################################

PACKAGE = vtkhdf
URL = https://github.com/johodges/fds.git
PACKAGEVERSION = git

MODE = mpi
MODULES = phdf5 

include ${MAKEINCLUDES}/Make.info

##
## cmake 
##
CMAKEFLAGS = -D CMAKE_BUILD_TYPE=Release \
    -D VTKHDF_USE_Fortran=OFF \
    -D VTKHDF_USE_AWSSDK=OFF
HASBIN = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

GITREPO = https://github.com/johodges/fds.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

