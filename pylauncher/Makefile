################################################################
####
#### Makefile for PYLAUNCHER installation
####
################################################################

PACKAGE = PYLAUNCHER
PACKAGEVERSION = 4.7
URL = https://github.com/TACC/pylauncher
MODE = core

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
NOLIB = 1
NOINC = 1
.PHONY: default_install
info ::
	@echo "make default_install : copy to install dir, make module"
default_install :
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}"  "${PACKAGEVERSION}"  "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}"  "${INSTALLROOT}"  "${INSTALLEXT}"  "${INSTALLVARIANT}" \
	    "${HOMEDIR}"  "${BUILDDIRROOT}"  "${MODE}" \
	    "${PREFIXOPTION}"  "${PREFIXEXTRA}" \
	 && mkdir -p $${prefixdir} \
	 && cp -r $${srcdir}/{docs,examples,README.md,scripts,src,tests} $${prefixdir}/ \
	 && make --no-print-directory varsmodule \
	      NOLIB=1 NOINC=1 \
	      DOCURL=https://docs.tacc.utexas.edu/software/pylauncher/ \
	      PYTHONPATHADD=src/pylauncher

include ${MAKEINCLUDES}/Make.vars

TGZURL = https://github.com/TACC/pylauncher/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/TACC/pylauncher.git
