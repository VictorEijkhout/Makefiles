################################################################
####
#### Makefile for Wannier90 installation
####
#### https://gitlab.com/wannier90-project/wannier90/-/issues/397
####
#### libfdf cloning in F jail: https://gitlab.com/wannier90-project/wannier90/-/issues/417
####
################################################################

PACKAGE = wannier90
URL = https://github.com/wannier-developers/wannier90
PACKAGEVERSION = 3.1.0

MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "make default_install : make-based install"
.PHONY: default_install
default_install :
	@scriptdir=$$(pwd) \
	 && source ${MAKEINCLUDES}/names.sh \
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
	     && rm -rf src/{obj,objp} \
	     && cp $${scriptdir}/make.inc.${TACC_FAMILY_COMPILER} ./make.inc \
	     && export DESTDIR="$${prefixdir}" && export PREFIX="" \
	         && make -j ${JCOUNT} all && make install \
	         && cp -r examples doc $${prefixdir}/ \
	         && mkdir -p $${prefixdir}/include \
	         && cp src/obj/*.mod $${prefixdir}/include \
	         && cp -r examples doc $${prefixdir}/ \
	    ) 2>&1 | tee $$installlog \
	 && make --no-print-directory varsmodule HASBIN=1 
include ${MAKEINCLUDES}/Make.vars
include ${MAKEINCLUDES}/Make.public


TGZURL = https://github.com/wannier-developers/wannier90/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://gitlab.com/wannier90-project/wannier90.git
SUBMODULE = 1
include ${MAKEINCLUDES}/Make.git
