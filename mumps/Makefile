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
MODULES = ${BASIC_MODULES} scotch metis

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
	    MODULENAME=pmumps MODULEALTNAME=mumps \
	    MODULES="${BASIC_MODULES} ptscotch parmetis"
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
	 && ( \
	     echo "Start build with PARMODE=${PARMODE}" \
	     && reportcompilers \
	     && requirenonzero builddir && rm -rf $$builddir && mkdir -p $$builddir \
	     && requirenonzero prefixdir && rm -rf $${prefixdir} && mkdir -p $${prefixdir} \
	     && cd $$srcdir \
	     && echo && echo "disable example making" && echo \
	     && sed -i Makefile -e '/cd examples/s/;/; echo/p' \
	     && grep examples Makefile \
	     && cp $${scriptdir}/Make.inc.${TACC_FAMILY_COMPILER}.${PARMODE} ./Makefile.inc \
	     && export DESTDIR="$${prefixdir}" && export PREFIX="" \
	         && make -j ${JCOUNT} \
	         && cp Makefile.inc    $${prefixdir}/ \
	         && cp -r include lib  $${prefixdir}/ \
	         && ( cd examples && make --no-print-directory clean ) \
	         && cp -r examples doc $${prefixdir}/ \
	         && cp -r SCILAB/examples $${prefixdir}/examples/SCILIB \
	    ) 2>&1 | tee $$installlog \
	 && make --no-print-directory varsmodule

include ${MAKEINCLUDES}/Make.vars
include ${MAKEINCLUDES}/Make.public

##
## download
##
TGZURL = https://mumps-solver.org/MUMPS_${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
