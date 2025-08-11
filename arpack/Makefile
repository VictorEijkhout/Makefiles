################################################################
####
#### Makefile for Arpack installation
#### 
################################################################

PACKAGE = ARPACK
ABOUT = Collection of Fortran77 subroutines designed to solve large scale eigenvalue problems.
URL = https://github.com/opencollab/arpack-ng
PACKAGEVERSION = 3.9.1
MODE = seq
MODULES = eigen mkl

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGSBASIC = \
    -D EXAMPLES=ON -D ICB=ON -D EIGEN=ON
ifeq "${TACC_SYSTEM}" "vista"
  CMAKEFLAGSBASIC += \
    -D LAPACK_LIBRARIES='${TACC_NVPL_LIB}/libnvpl_lapack_lp64_seq.so;${TACC_NVPL_LIB}/libnvpl_lapack_core.so' \
    -D BLAS_LIBRARIES='${TACC_NVPL_LIB}/libnvpl_blas_lp64_seq.so;${TACC_NVPL_LIB}/libnvpl_blas_core.so'
endif
CMAKEFLAGS = ${CMAKEFLAGSBASIC}
PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
CPTOINSTALLDIR = EXAMPLES

info ::
	@echo "make par : parallel installation"
	@echo "make seq : sequential installation"
ARPACKPARALLEL = OFF
.PHONY: seq par
seq :
	@make configure build JCOUNT=${JCOUNT} PACKAGEVERSION=${PACKAGEVERSION} \
	    CMAKEFLAGS="${CMAKEFLAGSBASIC}" \
	    ARPACKPARALLEL=OFF MODE=seq MODULENAME=arpack
par :
	@make configure build JCOUNT=${JCOUNT} PACKAGEVERSION=${PACKAGEVERSION} \
	    CMAKEFLAGS="${CMAKEFLAGSBASIC} -D MPI=ON" \
	    ARPACKPARALLEL=ON  MODE=mpi MODULENAME=parpack MODULENAMEALT=arpack
info ::
	@echo "make default_install : sequential version"
default_install :: seq 

##
## download
##
TGZURL = https://github.com/opencollab/arpack-ng/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

##
## git clone
##
GITREPO = https://github.com/opencollab/arpack-ng.git
include ${MAKEINCLUDES}/Make.git

##
## Clean 
##
include ${MAKEINCLUDES}/Make.clean
