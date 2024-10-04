################################################################
####
#### Makefile for Ncl installation
####
################################################################

PACKAGE = CESM2
PACKAGEVERSION = 2.1.3
URL = https://github.com/ESCOMP/CESM/tree/release-cesm2.2.2
MODE = mpi

include ${MAKEINCLUDES}/Make.info

.PHONY: default_install configure
default_install : configure 
configure : 
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 && cd $${srcdir} \
	 && git tag

HASBIN = 1
include ${MAKEINCLUDES}/Make.vars

TGZURL = https://github.com/ESCOMP/CESM/archive/refs/tags/release-cesm${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = 
include ${MAKEINCLUDES}/Make.git

