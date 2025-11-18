################################################################
####
#### Makefile for Eigen installation
####
#### eigen can not be used in Arpack with nvidia compiler:
#### https://gitlab.com/libeigen/eigen/-/issues/2846
####
################################################################

PACKAGE = EIGEN
ABOUT = C++ template library for linear algebra
PACKAGEVERSION = 5.0.0
URL = https://eigen.tuxfamily.org/index.php?title=Main_Page \
    https://gitlab.com/libeigen/eigen/-/tree/master
MODULES = 
MODE = seq

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
NOLIB = 1
PKGCONFIG = share/pkgconfig
CMAKEPREFIXPATHSET = 1
INCLUDELOC = include/eigen3
include ${MAKEINCLUDES}/Make.cbuild

.PHONY: default_install
default_install : configure build

TGZURL = https://gitlab.com/libeigen/eigen/-/archive/${PACKAGEVERSION}/eigen-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://gitlab.com/libeigen/eigen.git
include ${MAKEINCLUDES}/Make.git

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
