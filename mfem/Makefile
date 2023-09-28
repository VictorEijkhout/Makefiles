################################################################
####
#### Makefile for Mfem installation
####
#### only needed for the CHyPS CSA project
####
################################################################

PACKAGE = mfem
URL = https://github.com/mfem/mfem
PACKAGEVERSION = 4.4
## 4.5.2 has breaking changes for CHyPS.

MODE = mpi
MODULES = hypre64 metis

include ${MAKEINCLUDES}/Make.info

##
## cmake 
##
CMAKEFLAGS = -Wno-dev \
    -D CMAKE_BUILD_TYPE=Release \
    -D MFEM_USE_MPI=YES \
    -D MFEM_USE_LAPACK=YES \
    -D MFEM_DEBUG=NO \
    -D HYPRE_DIR=${TACC_HYPRE_DIR} \
    -D MFEM_USE_METIS_5=ON -D METIS_DIR=${TACC_METIS_DIR} \
    -D MFEM_USE_ADIOS2=ON -D ADIOS2_DIR=${TACC_ADIOS2_DIR}

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/mfem/mfem/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = git@github.com:mfem/mfem.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

