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

SCOTCHPARALLEL = OFF
CMAKEFLAGS = 
include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
CMAKEPREFIXPATHSET = 1
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
include ${MAKEINCLUDES}/Make.download

##
## git clone
##
GITREPO = https://gitlab.inria.fr/scotch/scotch.git
include ${MAKEINCLUDES}/Make.git

##
## clean
##
include ${MAKEINCLUDES}/Make.clean
