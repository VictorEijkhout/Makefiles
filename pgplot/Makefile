################################################################
####
#### Makefile for PGPLOT installation
####
################################################################

PACKAGE = PGPLOT
URL = https://sites.astro.caltech.edu/~tjp/pgplot/
PACKAGEVERSION = 5.2
MODE = seq
MODULES = 
include ${MAKEINCLUDES}/Make.info

#include ${MAKEINCLUDES}/Make.configure
#include ${MAKEINCLUDES}/Make.install

.PHONY: configure build
configure :: modules
	@echo "there is no configure rule"
build ::
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}"  "${PACKAGEVERSION}"  "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}"  "${INSTALLROOT}"  "${INSTALLEXT}"  "${INSTALLVARIANT}" \
	    "${HOMEDIR}"  "${BUILDDIRROOT}"  "${MODE}" \
	    "${PREFIXOPTION}"  "${PREFIXEXTRA}" \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 && requirenonzero srcdir \
	 && case "${TACC_FAMILY_COMPILER}" in \
	      ( "gcc" ) export conf=g90 ;; \
	      ( * ) echo "Only gcc implemented" && exit 1 ;; \
	    esac \
	 && echo " .. modifications in the source tree" \
	 && cp $${conf}.conf $${srcdir}/sys_linux/ \
	 && sed -e '/MFDRIV/s/!/ /' -e '/PNDRIV/s/!/ /' -e '/WDDRIV/s/!/ /' \
	        -i $${srcdir}/drivers.list \
	 && ( \
	    requirenonzero builddir && rm -rf $$builddir && mkdir -p $$builddir \
	     && requirenonzero prefixdir && mkdir -p $$prefixdir && rm -rf $$prefixdir/* \
	     && cd $${builddir} \
	     && $${srcdir}/makemake $${srcdir} linux $${conf} \
	     && make \
	    ) 2>&1 | tee $$installlog \

.PHONY: default_install
default_install : configure build 

TGZURL = ftp://ftp.astro.caltech.edu/pub/pgplot/pgplot${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean
