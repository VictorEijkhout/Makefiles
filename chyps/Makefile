################################################################
####
#### Makefile for Chyps installation
####
#### only needed for the CHyPS CSA project
####
################################################################

PACKAGE = chyps
URL = 
PACKAGEVERSION = git

MODE = mpi
MODULES = adios2 hypre metis mfem mutationpp precice petsc

include ${MAKEINCLUDES}/Make.info

##
## cmake 
##
CMAKEFLAGS = -DBUILD_SHARED_LIBS=ON \
    -D CMAKE_BUILD_TYPE=Release \
    -D BUILD_TESTING=OFF \
    -D Hypre_DIR=${TACC_HYPRE_DIR} \
    -D Metis_DIR=${TACC_METIS_DIR} \
    -D Mfem_DIR=${TACC_MFEM_DIR} \
    -D ADIOS2_DIR=${TACC_ADIOS2_DIR} \
    -D mutation++_DIR=${TACC_MUTATIONPP_DIR} \
    -D Eigen3_DIR=${TACC_EIGEN_DIR} \
    -D Precice_DIR=${TACC_PRECICE_DIR} \
    -D Petsc_DIR=${PETSC_DIR} \
    -D Petsc_LIBRARY=${TACC_PETSC_LIB}/libpetsc.so

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

GITREPO = git@github.com:chess-uiuc/chyps.git
SUBMODULE = 1
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

