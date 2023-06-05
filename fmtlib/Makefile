################################################################
####
#### Makefile for fmtlib
####
#### https://github.com/fmtlib
####
################################################################

PACKAGE = fmtlib
PACKAGEVERSION = 10.0.0
MODE = seq

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
ifeq "${TACC_FAMILY_COMPILER}" "intel"
  LIBDIR = lib64
  PKGCONFIG = lib64/pkgconfig
else
  PKGCONFIG = lib/pkgconfig
endif
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make docs, docview (after configure)"
.PHONY: docs docview
docs :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setnames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && if [ ! -d "$${builddir}" ] ; then \
	      echo "Please rerun configuration: need build dir" && exit 1 ; fi \
	 && cd "$${builddir}" && make doc
docview :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setnames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && cd $${installdir}/doc/html && open index.html

TGZURL = https://github.com/fmtlib/fmt/archive/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
