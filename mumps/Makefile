################################################################
####
#### Makefile for MUMPS installation
####
################################################################

PACKAGE = MUMPS
PACKAGEVERSION = 5.8.0
URL = https://www.mumps.org/
MODE = mpi
BASIC_MODULES = mkl
MODULES = ${BASIC_MODULES} scotch

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "make par, seq : makefile-based install"
	@echo "make default_install : par"
.PHONY: seq par default_install
PARMODE = par
seq :
	@${MAKE} default_install MODE=seq PARMODE=seq
par :
	@${MAKE} default_install MODE=mpi PARMODE=par \
	    MODULES="${BASIC_MODULES} ptscotch"
default_install :
	@scriptdir=$$(pwd) \
	 && source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 && if [ "${PARMODE}" = "seq" ] ; then setcompilers ; else setmpicompilers ; fi \
	 && reportcompilers \
	 && ( \
	    requirenonzero builddir && rm -rf $$builddir && mkdir -p $$builddir \
	     && requirenonzero prefixdir && rm -rf $${prefixdir} && mkdir -p $${prefixdir}/lib \
	     && cd $$srcdir \
	     && cp $${scriptdir}/Make.inc.${TACC_FAMILY_COMPILER}.${PARMODE} ./Makefile.inc \
	     && export DESTDIR="$${prefixdir}" && export PREFIX="" \
	         && make -j ${JCOUNT} \
	         && cp -r include lib  $${prefixdir}/ \
	         && cp -r examples doc $${prefixdir}/ \
	    ) 2>&1 | tee $$installlog \
	 && make --no-print-directory varsmodule HASBIN=1 
include ${MAKEINCLUDES}/Make.vars
include ${MAKEINCLUDES}/Make.public

##
## download
##
TGZURL = https://mumps-solver.org/MUMPS_${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
