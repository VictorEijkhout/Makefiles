################################################################
####
#### Makefile for fmtlib
####
################################################################

PACKAGE = fmtlib
URL = https://github.com/fmtlib/fmt
PACKAGEVERSION = 11.0.2
MODE = seq

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
CMAKEFLAGS = \
    -D FMT_TEST=OFF
PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
include ${MAKEINCLUDES}/Make.cbuild
default_install : configure build

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
