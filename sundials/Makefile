################################################################
####
#### Makefile for Sundials installation
####
################################################################

PACKAGE = SUNDIALS
ABOUT = SUite of Nonlinear and DIfferential/ALgebraic equation Solvers.
URL = https://computing.llnl.gov/projects/sundials \
    https://github.com/LLNL/sundials \
    https://sundials.readthedocs.io/en/latest/sundials/Install_link.html#installation-cmake
PACKAGEVERSION = 7.4.0
MODE = mpi

include ${MAKEINCLUDES}/Make.info

PARALLEL = ON
CMAKEFLAGS = \
    -D ENABLE_MPI=${PARALLEL} \
    -D BUILD_FORTRAN_MODULE_INTERFACE=ON
CMAKEPREFIXPATHSET = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
info ::
	@echo "make default_install : plain configure build"
default_install : par
par: configure build
seq :
	@make --no-print-directory configure build \
	    JCOUNT=${JCOUNT} PACKAGEVERSION=${PACKAGEVERSION} \
	    MODE=seq PARALLEL=OFF \
	    MODULENAME=ssundials MODULENAMEALT=sundials

HASBIN = 1

TGZURL = https://github.com/LLNL/sundials/releases/download/v${PACKAGEVERSION}/sundials-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean
