################################################################
####
#### Makefile for Ncl installation
####
################################################################

PACKAGE = ncl
PACKAGEVERSION = 6.6.2
PACKAGEVERSION_NODOT = 662
URL = https://www.ncl.ucar.edu/ https://www.ncl.ucar.edu/Download/build_from_src.shtml
MODE = seq

include ${MAKEINCLUDES}/Make.info

.PHONY: mkinfo
mkinfo :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 && cd $$srcdir && make Info


CONFIGUREFLAGS = 
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install

TGZURL = https://www.earthsystemgrid.org/dataset/ncl.${PACKAGEVERSION_NODOT}.src/file/ncl_ncarg-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/flame/ncl.git
include ${MAKEINCLUDES}/Make.git

