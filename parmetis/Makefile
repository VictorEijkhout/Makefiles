################################################################
####
#### Makefile for Parmetis installation
####
#### this is the 4.0.3 version from
#### https://github.com/scibuilder/parmetis
#### this builds gklib & metis too 
####
################################################################

PACKAGE = parmetis
PACKAGEVERSION = 4.0.3
MODE = mpi
MODULES = 
#gklib metis

include ${MAKEINCLUDES}/Make.info
CMAKEFLAGS=\
    -DCMAKE_C_FLAGS=\"-O2 -qopenmp\" 

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

include ${MAKEINCLUDES}/Make.clean

GITREPO = https://github.com/scibuilder/parmetis.git
include ${MAKEINCLUDES}/Make.git
ZIPURL = https://github.com/scibuilder/parmetis/archive/refs/heads/master.zip
include ${MAKEINCLUDES}/Make.download
