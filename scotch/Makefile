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
PACKAGEVERSION = 7.0.7
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

SCOTCHPARALLEL = OFF
INTSIZE = 32
CMAKEFLAGS = \
    -D BUILD_PTSCOTCH:BOOL=${SCOTCHPARALLEL} \
    -D INTSIZE:STRING=${INTSIZE} \
    -D BUILD_FORTRAN=ON
include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
CMAKEPREFIXPATHSET = 1
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "================ scotch rules"
	@echo "make par : parallel installation"
	@echo "make seq : sequential installation"
.PHONY: seq32 par32
seq32 :
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    PACKAGEVERSION=${PACKAGEVERSION} \
	    SCOTCHPARALLEL=OFF MODE=seq \
	    INTSIZE=32 INSTALLVARIANT=i32 INSTALLEXT=i32 \
	    MODULENAME=scotch
par32 :
	@if [ "${TACC_FAMILY_COMPILER}" = "gcc" ] ; then \
	    export I_MPI_FCFLAGS="-I${TACC_IMPI_INC}/mpi/gfortran/11.1.0" \
	 ; fi \
	 && make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    PACKAGEVERSION=${PACKAGEVERSION} \
	    SCOTCHPARALLEL=ON  MODE=mpi \
	    INTSIZE=32 INSTALLVARIANT=i32 INSTALLEXT=i32 \
	    MODULENAME=ptscotch 
.PHONY: seq64 par64
seq64 :
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    PACKAGEVERSION=${PACKAGEVERSION} \
	    SCOTCHPARALLEL=OFF MODE=seq \
	    INTSIZE=64 INSTALLVARIANT=i64 INSTALLEXT=i64 \
	    MODULENAME=scotch
par64 :
	@if [ "${TACC_FAMILY_COMPILER}" = "gcc" ] ; then \
	    export I_MPI_FCFLAGS="-I${TACC_IMPI_INC}/mpi/gfortran/11.1.0" \
	 ; fi \
	 && make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    PACKAGEVERSION=${PACKAGEVERSION} \
	    SCOTCHPARALLEL=ON  MODE=mpi \
	    INTSIZE=64 INSTALLVARIANT=i64 INSTALLEXT=i64 \
	    MODULENAME=ptscotch 
info ::
	@echo "make default_install : sequential and parallel versions"
default_install :: 
	@make --no-print-directory JCOUNT=${JCOUNT} PACKAGEVERSION=${PACKAGEVERSION} \
	     seq32 par32

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
