################################################################
####
#### Makefile for Ncl installation
####
################################################################

PACKAGE = CESM2
ABOUT = The Community Earth System Model
PACKAGEVERSION = 2.1.3
URL = https://github.com/ESCOMP/CESM/tree/release-cesm2.2.2
MODE = mpi

include ${MAKEINCLUDES}/Make.info

.PHONY: default_install configure
default_install : configure 
configure : 
	@source ${MAKEINCLUDES}/names.sh \
	 && PACKAGEVERSION=git \
	    setnames \
	    "${PACKAGE}"  "${PACKAGEVERSION}"  "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}"  "${INSTALLROOT}"  "${INSTALLEXT}"  "${INSTALLVARIANT}" \
	    "${HOMEDIR}"  "${BUILDDIRROOT}"  "${MODE}" \
	    "${PREFIXOPTION}"  "${PREFIXEXTRA}" \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 && cd $${srcdir} \
	 && git tag

HASBIN = 1
include ${MAKEINCLUDES}/Make.vars

DOWNLOADURL = https://github.com/ESCOMP/CESM/archive/refs/tags/release-cesm${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/escomp/cesm.git
include ${MAKEINCLUDES}/Make.git

