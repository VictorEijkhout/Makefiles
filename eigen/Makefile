################################################################
####
#### Makefile for Eigen installation
####
################################################################

PACKAGE = EIGEN
PACKAGEVERSION = 3.4.0
URL = https://eigen.tuxfamily.org/index.php?title=Main_Page
MODULES = 
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIG = share/pkgconfig
CMAKEPREFIXPATHSET = 1
include ${MAKEINCLUDES}/Make.cmake
NOLIB = 1
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install : configure build

TGZURL = https://gitlab.com/libeigen/eigen/-/archive/${PACKAGEVERSION}/eigen-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

.PHONY: docs
docs :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 && if [ ! -d $$builddir ] ; then \
	        echo "Not finding build dir: $$builddir" && exit 1 ; fi \
	 && cd $$builddir/doc \
	 && doxygen Doxyfile \
	 && cp -r html $$homedir/$${package}-$${packageversion}-manual
not_latex_docs :
	foo \
	 && ( cat Doxyfile ; echo "GENERATE_LATEX=YES" ) | doxygen - \
	 && cd latex && pdflatex refman && pdflatex refman \
	 && cp refman.pdf $$homedir/$${package}-$${packageversion}-refman.pdf
