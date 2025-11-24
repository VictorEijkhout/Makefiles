################################################################
####
#### Makefile for mpl
####
################################################################

PACKAGE = MPL
ABOUT = A C++17 message passing library based on MPI
PACKAGEVERSION = 0.3.0
URL = https://github.com/rabauke/mpl
MODE = mpi

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install : configure build
CMAKEPREFIXPATHSET = 1
MODULE_HELP = 1
module_help :
	@for line in \
	    "Usage:" \
	    "  import mpl" \
	    "and use" \
	    "   #include <mpl/mpl.hpp>" \
	    "in your code" \
	; do echo "$$line" ; done 


info ::
	@echo "================ docs rules"
	@echo "make docs : doxygen html & latex docs"
	@echo "    [ LATEXDOC=0/1 (default=${LATEXDOC}) ]"
.PHONY: docs opendocs
LATEXDOC=0
docs :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
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
info ::
	@echo "make opendocs"
opendocs :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 && doctarget=$${installdir}/../docs-${PACKAGEVERSION} \
	 && ls $$doctarget \
	 && open $$doctarget/html/index.html
info ::
	@echo "make findexample ROUTINE=...."
	@echo "    [ FLAGS= ]"
ROUTINE = allreduce
.PHONY: findexample
findexample :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setdirlognames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} ${INSTALLEXT} \
	 && grep ${FLAGS} "${ROUTINE}" $$srcdir/examples/*.cc

GITREPO = git@github.com:rabauke/mpl.git
include ${MAKEINCLUDES}/Make.git

DOWNLOADURL = https://github.com/rabauke/mpl/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
