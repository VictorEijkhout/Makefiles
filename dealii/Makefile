################################################################
####
#### Makefile for Dealii installation
####
#### Kokkos problem:
#### https://github.com/dealii/dealii/issues/18928
####
################################################################

PACKAGE = DEALII
ABOUT = Finite Element Library
URL = https://github.com/dealii/dealii
PACKAGEVERSION = 9.7.1
MODE = mpi
MODULES = p4est petsc phdf5 sundials trilinos \
    blaslapack
# intel 23 problem solved in 24
# support ticket 05936427

CONFIGURATION = default
include ${MAKEINCLUDES}/Make.info

##
## https://dealii.org/developer/readme.html
##
COMPLEX = OFF
CMAKEFLAGS = \
    -D DEAL_II_WITH_MPI=ON \
    -D CMAKE_CXX_FLAGS="-fopenmp" \
    -D DEAL_II_WITH_LAPACK=ON \
    -D DEAL_II_WITH_TASKFLOW=OFF \
    -D DEAL_II_COMPONENT_EXAMPLES=OFF \
    -D DEAL_II_WITH_PETSC_COMPLEX=${COMPLEX} \
    -D DEAL_II_WITH_COMPLEX_VALUES=${COMPLEX} \
    -D DEAL_II_WITH_TBB=OFF \
    -D DEAL_II_WITH_HDF5=ON \
    -D HDF5_DIR=${TACC_PHDF5_DIR} \
    -D DEAL_II_WITH_P4EST=ON    -D P4EST_DIR=${TACC_P4EST_DIR} \
    -D DEAL_II_WITH_PETSC=ON    -D PETSC_DIR=${PETSC_DIR} -D PETSC_ARCH=${PETSC_ARCH} \
    -D DEAL_II_WITH_SUNDIALS=ON -D SUNDIALS_DIR=${TACC_SUNDIALS_DIR} \
    -D DEAL_II_WITH_TRILINOS=ON -D TRILINOS_DIR=${TACC_TRILINOS_DIR}

##
## Lapack
## https://dealii.org/developer/external-libs/lapack.html
##
ifdef TACC_MKL_DIR
  CMAKEFLAGS += \
    -D LAPACK_DIR=${TACC_MKL_LIB}
endif

# ifeq "${TACC_FAMILY_COMPILER}" "intel"
# CMAKEFLAGS += \
#     -D DEAL_II_CXX_FLAGS='-Wno-\#warnings -save-temps -v' \
#     -D DEAL_II_LINKER_FLAGS=-fuse-ld=bfd 
# endif

## was this always wrong or did it go away? PKGCONFIG = lib/pkgconfig
CMAKEPREFIXPATHSET = 1
CPTOINSTALLDIR = examples
# https://github.com/dealii/dealii/issues/18254
CMAKEBUILDTYPE = DebugRelease
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
# https://github.com/dealii/dealii/issues/18254
real: 
	@cd ${TACC_PETSC_DIR} \
	 && has=$$( grep "define PETSC_USE_COMPLEX" include/petscconf.h | cut -d ' ' -f 3 ) \
	 && if [ ! -z "$${has}" -a "$$has" -ne 0 ] ; then \
	      echo "Real version needs real PETSC" && exit 1 ; fi
	@make --no-print-directory configure build PACKAGEVERSION=${PACKAGEVERSION} \
	    INSTALLEXT="real" INSTALLVARIANT="real" COMPLEX=OFF
complex :
	@cd ${TACC_PETSC_DIR} \
	 && has=$$( grep "define PETSC_USE_COMPLEX" include/petscconf.h | cut -d ' ' -f 3 ) \
	 && if [ -z "$${has}" -o "$$has" = "0" ] ; then \
	      echo "Complex version needs complex PETSC" && exit 1 ; fi
	@make --no-print-directory configure build PACKAGEVERSION=${PACKAGEVERSION} \
	    INSTALLEXT="complex" INSTALLVARIANT="complex" COMPLEX=ON

.PHONY: default_install real complex
info ::
	@echo "make default_install : real version, there is also complex"
# https://github.com/dealii/dealii/issues/18976
default_install : real

# need modules loaded for shared libs
DEPENDSON = p4est sundials trilinos
# needs to find the hdf5 libs
DEPENDSONCURRENT = phdf5

GITREPO = https://github.com/dealii/dealii/archive/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.git

DOWNLOADURL = https://github.com/dealii/dealii/releases/download/v${PACKAGEVERSION}/dealii-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

##
## testing
##
info ::
	@echo && echo "================" && echo "make test [ TESTSTEP=nnn ]"
.PHONY: test
TESTSTEP = 17
test :
	@cd testing-${PACKAGEVERSION} \
	 && make --no-print-directory test STEP=${TESTSTEP}

