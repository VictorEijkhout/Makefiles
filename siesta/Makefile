################################################################
####
#### Makefile for Siesta installation
####
#### Fortran debug flags: https://gitlab.com/siesta-project/siesta/-/issues/397
#### libfdf cloning in F jail: https://gitlab.com/siesta-project/siesta/-/issues/417
#### scalapack fortran interface: https://gitlab.com/siesta-project/siesta/-/issues/438
####
#### Vista:
#### https://gitlab.com/siesta-project/siesta/-/issues/441
####
################################################################

PACKAGE = siesta
URL = https://siesta-project.org/siesta/index.html
PACKAGEVERSION = 5.2.1

MODE = mpi
MODULES = fftw3 mkl 

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig

##
## cmake version
##
# https://docs.siesta-project.org/projects/siesta/en/stable/installation/build-issues.html#siesta-does-not-recognize-intel-mkl-libraries
## CMAKEPREP = find . -name \*make -exec sed -i -e 's/\-xHost//' -e 's/\-check//' {} \;
## CMAKEBUILDDEBUG = 1
CMAKEFLAGS = \
    -D SIESTA_WITH_FLOOK=FALSE \
    -D Fortran_FLAGS_DEBUG=-O0 \
    -D CMAKE_Fortran_FLAGS=-fPIE \
    \
    -D SIESTA_WITH_WANNIER90=ON \
    \
    -D SIESTA_WITH_FFTW=TRUE
ifdef TACC_MKL_DIR
  CMAKEFLAGS += \
    -D SCALAPACK_LIBRARY_DIR=${TACC_MKL_LIB} \
    -D SCALAPACK_LIBRARY='-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core'
else ifdef TACC_NVPL_DIR
  CMAKEFLAGS += \
    -D BLAS_LIBRARY_DIR=${TACC_NVPL_LIB} \
    -D BLAS_LIBRARY='-lnvpl_blas_lp64_seq -lnvpl_blas_core' \
    -D LAPACK_LIBRARY_DIR=${TACC_NVPL_LIB} \
    -D LAPACK_LIBRARY='-lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core' \
    -D SCALAPACK_LIBRARY_DIR=${TACC_NVPL_LIB} \
    -D SCALAPACK_LIBRARY='-lnvpl_scalapack_ilp64 -lnvpl_blacs_lp64_openmpi5'
endif

HASBIN = 1
CPTOINSTALLDIR = Examples
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build


TGZURL = https://gitlab.com/siesta-project/siesta/-/archive/${PACKAGEVERSION}/siesta-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://gitlab.com/siesta-project/siesta.git
SUBMODULE = 1
include ${MAKEINCLUDES}/Make.git
