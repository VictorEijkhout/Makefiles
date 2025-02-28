################################################################
####
#### Makefile for Scotch installation
####
#### note: also comes in C++ and Fortran download
#### 
#### see release_docs/INSTALL_CMake.txt
#### 
################################################################

PACKAGE = SCOTCH
URL = https://www.labri.fr/perso/pelegrin/scotch/ \
    https://gitlab.inria.fr/scotch/scotch
PACKAGEVERSION = 7.0.6
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

SCOTCHPARALLEL = ON
HDFFORTRAN = ON

ifeq "${TACC_SYSTEM}" "macbookair"
  LIBEXT = dylib
else
  LIBEXT = so
endif

#### test failure with NVidia compiler:
#### https://hdfgroup.atlassian.net/servicedesk/customer/portal/6/HELP-2539
TESTING = ON

CMAKEFLAGS = 
include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
FAMILY = hf5
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "================ scotch rules"
	@echo "make par : parallel installation"
	@echo "make seq : sequential installation"
.PHONY: seq par
seq :
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    PACKAGEVERSION=${PACKAGEVERSION} \
	    SCOTCHPARALLEL=OFF MODE=seq \
	    MODULENAME=scotch
par :
	@if [ "${TACC_FAMILY_COMPILER}" = "gcc" ] ; then \
	    export I_MPI_FCFLAGS="-I${TACC_IMPI_INC}/mpi/gfortran/11.1.0" \
	 ; fi \
	 && make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    PACKAGEVERSION=${PACKAGEVERSION} \
	    SCOTCHPARALLEL=ON  MODE=mpi \
	    MODULENAME=ptscotch 
info ::
	@echo "make default_install : sequential and parallel versions"
default_install :: 
	@make --no-print-directory seq par JCOUNT=${JCOUNT} PACKAGEVERSION=${PACKAGEVERSION} \
	    seq par

##
## download
##
TGZURL = https://gitlab.inria.fr/scotch/scotch/-/archive/v${PACKAGEVERSION}/scotch-v${PACKAGEVERSION}.tar.gz

##
## git clone
##
GITREPO = https://gitlab.inria.fr/scotch/scotch.git
include ${MAKEINCLUDES}/Make.git

##
## clean
##
include ${MAKEINCLUDES}/Make.clean
