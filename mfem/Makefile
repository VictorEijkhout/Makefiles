################################################################
####
#### Makefile for Mfem installation
####
#### only needed for the CHyPS CSA project
####
################################################################

PACKAGE = mfem
URL = https://github.com/mfem/mfem
PACKAGEVERSION = 4.7
## 4.5.2 has breaking changes for CHyPS.
HYPREVERSION = 2.31.0

MODE = mpi
MODULES = adios2 hypre/${HYPREVERSION}-i64 metis

include ${MAKEINCLUDES}/Make.info

##
## cmake 
##
CMAKEFLAGS = -Wno-dev \
    -D CMAKE_BUILD_TYPE=Release \
    -D MFEM_USE_MPI=YES \
    -D MFEM_USE_LAPACK=YES \
    -D LAPACK_LIBRARIES=${TACC_NVPL_LIB}/libnvpl_lapack_lp64_seq.so \
    -D BLAS_LIBRARIES=${TACC_NVPL_LIB}/libnvpl_blas_lp64_seq.so \
    -D MFEM_DEBUG=NO \
    -D HYPRE_DIR=${TACC_HYPRE_DIR} \
    -D MFEM_USE_METIS_5=ON -D METIS_DIR=${TACC_METIS_DIR} \
    -D MFEM_USE_ADIOS2=ON -D ADIOS2_DIR=${TACC_ADIOS2_DIR}

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

TGZURL = https://github.com/mfem/mfem/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = git@github.com:mfem/mfem.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean
