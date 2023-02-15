################################################################
####
#### Makefile for mpl
####
################################################################

PACKAGE = MPL
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 0.2.0
MODE = mpi

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
# EXTRABUILDTARGETS=documentation
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make docs : doxygen html & latex docs"
	@echo "    [ LATEXDOC=0/1 (default=${LATEXDOC}) ]"
	@echo "make opendocs"
.PHONY: docs opendocs
LATEXDOC=0
docs :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} ${INSTALLEXT} \
	 && cd $${builddir}/doc/doxygen \
	 && doctarget=$${installdir}/../docs-${PACKAGEVERSION} \
	 && if "${LATEXDOC}" = "1" ] ; then \
	        LATEXYESNO=YES ; else LATEXYESNO=NO ; fi \
	 && ( cat Doxyfile && echo GENERATE_LATEX=${LATEXYESNO} ) | doxygen - \
	 && if [ "${LATEXDOC}" = "1" ] ; then \
	      ( cd latex \
	        && yes | tr y q | pdflatex refman || echo "skip errors" \
	        && yes | tr y q | pdflatex refman || echo "skip errors" \
	      ) ; fi \
	 && mkdir -p $$doctarget \
	 && for t in html latex ; do \
	      if [ -d $$t ] ; then \
	        if [ -d $${doctarget}/$$t ] ; then \
	          rm -rf $${doctarget}/$$t ; fi \
	         && cp -r $$t $${doctarget}/ \
	      ; fi \
	    ; done
opendocs :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} ${INSTALLEXT} \
	 && doctarget=$${installdir}/../docs-${PACKAGEVERSION} \
	 && ls $$doctarget \
	 && open $$doctarget/html/index.html

GITREPO = git@github.com:rabauke/mpl.git
include ${MAKEINCLUDES}/Make.git

TGZURL = https://github.com/rabauke/mpl/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
