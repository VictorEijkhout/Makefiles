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

# intel compiler considers missing proto as fatal
CMAKEFLAGS = --compile-no-warning-as-error
include ${MAKEINCLUDES}/Make.cmake
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
