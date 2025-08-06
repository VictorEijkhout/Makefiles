################################################################
####
#### Makefile for Rmt installation
####
################################################################

PACKAGE = rmt
ABOUT = RMT (R-matrix with time-dependence) is part of the UK-AMOR suite of codes for tackling problems in Atomic Molecular and Optical physics problems with the Rmatrix approach.
URL = https://gitlab.com/Uk-amor/RMT/rmt
PACKAGEVERSION = git

MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

##
## cmake version
##
include ${MAKEINCLUDES}/Make.cmake
CMAKESUBDIR = source
NOINSTALL = 1
include ${MAKEINCLUDES}/Make.cbuild
build ::
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 && mkdir -p $$installdir \
	 && cp -r $${builddir}/{bin,lib,modules,programs} $${installdir}/ \
	 && chmod -R g+rX,o+rX $${installdir}
.PHONY: default_install
default_install : configure build

GITREPO = https://gitlab.com/Uk-amor/RMT/rmt.git
include ${MAKEINCLUDES}/Make.git
