################################################################
####
#### Makefile for SCONS installation
####
################################################################

PACKAGE = SCONS
ABOUT = Powerful Software Build System
PACKAGEVERSION = 4.9.1
URL = https://scons.org/ \
    https://github.com/SCons/scons
DOCURL = 
MODE = core

include ${MAKEINCLUDES}/Make.info

NOLIB = 1
NOINC = 1
HASBIN = 1
# include ${MAKEINCLUDES}/Make.scons
# include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : build varsmodule public"
.PHONY: build default_install
default_install : build varsmodule public
build ::
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && requirenonzero prefixdir && mkdir -p $${prefixdir} && rm -rf $${prefixdir}/* \
	 && requirenonzero builddir && rm -rf $$builddir && mkdir -p $$builddir \
	 && requirenonzero logfilesdir \
	 && echo "logfiles dir: $${logfilesdir}." \
	 && mkdir -p $${logfilesdir} \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 \
	 && echo "using python3=$$( which python3 )" \
	 && python3 -m pip install --target "$${prefixdir}" \
	      scons
PYTHONPATHADDREL = 1
include ${MAKEINCLUDES}/Make.vars

DOWNLOADURL = https://github.com/SCons/scons/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
