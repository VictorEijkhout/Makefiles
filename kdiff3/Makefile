################################################################
####
#### Makefile for KDIFF3 installation
####
################################################################

PACKAGE = KDIFF3
ABOUT = KDE file diff utility
URL = https://kdiff3.sourceforge.net/
PACKAGEVERSION = 0.9.98
MODE = core
MODULES = 

include ${MAKEINCLUDES}/Make.info

.PHONY: default_install configure build
configure ::
	@echo "there is no configure rule; go straight to default_install"
build ::
	@echo "there is no build rule; go straight to default_install"
default_install :
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
	     && requirenonzero prefixdir && rm -rf $$prefixdir/ \
	     && mkdir -p $$prefixdir/{bin,lib} \
	     && requirenonzero srcdir && cd $$srcdir \
	     && ./configure -o $${builddir} kde4 \
	    ) 2>&1 | tee $$installlog \
	 && make --no-print-directory varsmodule \
	      CMAKEPREFIXPATHSET=1
NOLIB = 1
HASBIN = 1
include ${MAKEINCLUDES}/Make.vars

TGZURL = https://sourceforge.net/projects/kdiff3/files/kdiff3/${PACKAGEVERSION}/kdiff3-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean
