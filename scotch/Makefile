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
# 7.0.8 is out
# Packagers are  encouraged to configure SCOTCH  so as to make the  widest use
# of threads, that is: -DCOMMON_PTHREAD -DSCOTCH_PTHREAD -DSCOTCH_PTHREAD_MPI,
# and also -DCOMMON_PTHREAD_AFFINITY_LINUX whenever possible.  On Windows-like
# environments, the -DCOMMON_THREAD_WIN32 flag is now available.
# The  -DSCOTCH_DETERMINISTIC  flag  may be  set or not,  depending on  former
# practice, given that users can now set this behavior at run time.
# Please refer to the file INSTALL.txt in the root directory.
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
	@echo "make seq32/64 : sequential version"
	@echo "make par32/64 : sequential version"
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
