################################################################
####
#### Makefile for Dealii installation
####
################################################################

PACKAGE = DEALII
ABOUT = Finite Element Library
URL = https://github.com/dealii/dealii
PACKAGEVERSION = 9.6.2
MODE = mpi
MODULES = p4est petsc sundials trilinos \
    mkl
# intel 23 problem solved in 24
# support ticket 05936427

CONFIGURATION = default
include ${MAKEINCLUDES}/Make.info

##
## https://dealii.org/developer/readme.html
##
CMAKEFLAGS = \
    -D DEAL_II_WITH_MPI=ON \
    -D DEAL_II_WITH_LAPACK=ON \
    -D DEAL_II_COMPONENT_EXAMPLES=OFF \
    -D DEAL_II_WITH_TBB=OFF \
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
# https://github.com/dealii/dealii/issues/18254
CMAKEBUILDTYPE = DebugRelease
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: real complex
info ::
	@echo "================"
	@echo "make real : default setup"
	@echo "make complex: configure build with complex. Make sure to load complex petsc!"
real :
	@make --no-print-directory configure build PACKAGEVERSION=${PACKAGEVERSION}
complex :
	@make --no-print-directory configure build PACKAGEVERSION=${PACKAGEVERSION} \
	    INSTALLEXT="complex" CMAKEFLAGS="${CMAKEFLAGS} -DDEAL_II_WITH_PETSC_COMPLEX=ON"

.PHONY: default_install
info ::
	@echo "make default_install : real version, there is also complex"
default_install : real
# need to make sure we load real/complex petsc
DEPENDSONCURRENT = petsc

GITREPO = https://github.com/dealii/dealii/archive/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.git

TGZURL = https://github.com/dealii/dealii/releases/download/v${PACKAGEVERSION}/dealii-${PACKAGEVERSION}.tar.gz
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

