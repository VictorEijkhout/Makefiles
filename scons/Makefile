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
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : 
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
	 && python -m pip install --user --target "$${prefixdir}" \
	      scons

TGZURL = https://github.com/SCons/scons/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
