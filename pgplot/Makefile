################################################################
####
#### Makefile for PGPLOT installation
####
#### https://consult.tacc.utexas.edu/Ticket/Display.html?id=101229
####
################################################################

PACKAGE = PGPLOT
URL = https://sites.astro.caltech.edu/~tjp/pgplot/
PACKAGEVERSION = 5.2
MODE = seq
MODULES = mkl
include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.vars
include ${MAKEINCLUDES}/Make.public
#include ${MAKEINCLUDES}/Make.configure
#include ${MAKEINCLUDES}/Make.install
.PHONY: configure build
configure :: modules
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
	      ( "intel" ) export conf=intel ;; \
	      ( * ) echo "Not implemented for compiler=${TACC_FAMILY_COMPILER}" && exit 1 ;; \
	    esac \
	 && echo "Modifications in the source tree" \
	 && cp $${conf}.conf $${srcdir}/sys_linux/ \
	 && echo " .. skip PNDRIV MFDRIV for compilation problem" \
	 && cp drivers.list $${srcdir} \
	 && find $${srcdir} -name \*.f \
	      -exec sed -i \
	                -e "/'X/s/'\([0-9A-F]\+\)'X/int(Z'\1',KIND=KIND(2))/" \
	                {} \;

build ::
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}"  "${PACKAGEVERSION}"  "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}"  "${INSTALLROOT}"  "${INSTALLEXT}"  "${INSTALLVARIANT}" \
	    "${HOMEDIR}"  "${BUILDDIRROOT}"  "${MODE}" \
	    "${PREFIXOPTION}"  "${PREFIXEXTRA}" \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 && case "${TACC_FAMILY_COMPILER}" in \
	      ( "gcc" ) export conf=g90 ;; \
	      ( "intel" ) export conf=intel ;; \
	      ( * ) echo "Not implemented for compiler=${TACC_FAMILY_COMPILER}" && exit 1 ;; \
	    esac \
	 && ( \
	    builddir=$${srcdir}/build && rm -rf $$builddir && mkdir -p $$builddir \
	     && requirenonzero prefixdir && mkdir -p $$prefixdir && rm -rf $$prefixdir/* \
	     && echo " .. building in the build dir does not work, using src" \
	     && cd $${srcdir} \
	     && echo " .. make make" \
	     && $${srcdir}/makemake $${srcdir} linux $${conf} \
	     && echo " .. do make" \
	     && make \
	     && echo " .. installing by explicit copy" \
	     && mkdir -p $${prefixdir}/lib && cp *.{a,so} $${prefixdir}/lib \
	     && mkdir -p $${prefixdir}/share && cp grfont.dat rgb.txt $${prefixdir}/share \
	     && cp -r $${srcdir}/examples $${prefixdir} \
	    ) 2>&1 | tee $$installlog \
	 && make --no-print-directory varsmodule

.PHONY: default_install
default_install : configure build 

TGZURL = ftp://ftp.astro.caltech.edu/pub/pgplot/pgplot${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean

#   445 |       COMBUF(1) = int('A701'X,KIND=KIND(2))
#       |                  1
# Error: Arithmetic overflow converting INTEGER(4) to INTEGER(2) at (1). This check can be disabled with the option ‘-fno-range-check’
# /work2/00434/eijkhout/pgplot/pgplot-5.2/drivers/mfdriv.f:197:72:
