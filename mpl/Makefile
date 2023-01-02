################################################################
####
#### Makefile for mpl
####
################################################################

PACKAGE = MPL
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 0.2.0
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODE = seq

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
# EXTRABUILDTARGETS=documentation
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make docs : doxygen html & latex docs"
.PHONY: docs
LATEXDOC=NO
docs :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} ${INSTALLEXT} \
	 && cd $${builddir}/doc/doxygen \
	 && doctarget=$${installdir}/../docs-${PACKAGEVERSION} \
	 && mkdir -p $$doctarget \
	 && ( cat Doxyfile && echo GENERATE_LATEX=${LATEXDOC} ) | doxygen - \
	 && if [ "${LATEXDOC}" = "YES" ] ; then \
	      ( cd latex && pdflatex refman && pdflatex refman  ) ; fi \
	 && for t in html latex ; do \
	      if [ -d $${doctarget}/$$t ] ; then \
	        rm -rf $${doctarget}/$$t ; fi \
	       && cp -r $$t $${doctarget}/ \
	    ; done

GITREPO = git@github.com:rabauke/mpl.git
include ${MAKEINCLUDES}/Make.git

TGZURL = https://github.com/rabauke/mpl/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
