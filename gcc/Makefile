################################################################
####
#### Makefile for Gcc installation
####
################################################################

PACKAGE = gcc
PACKAGEVERSION = 14.2
URL = 
MODE = core
MODULES = 

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "================ bootstrap install"
	@echo "make build"
.PHONY: configure build
configure :: modules
	@echo "there is no configure rule"
TOOLSET = ${TACC_FAMILY_COMPILER}
build ::
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 && ( \
	    requirenonzero builddir && rm -rf $$builddir && mkdir -p $$builddir \
	     && requirenonzero prefixdir && mkdir -p $$prefixdir && rm -rf $$prefixdir/* \
	     && cd $$srcdir \
	    ) 2>&1 | tee $$installlog \
	 && make --no-print-directory varsmodule \
	      CMAKEPREFIXPATHSET=1
include ${MAKEINCLUDES}/Make.vars
## NOLIB = 1
include ${MAKEINCLUDES}/Make.public

GITREPO = https://gcc.gnu.org/git/gcc.git
BRANCH = releases/gcc-14.2
include ${MAKEINCLUDES}/Make.git
