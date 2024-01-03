################################################################
####
#### Makefile for GKlib installation
####
################################################################

PACKAGE = gklib
URL = https://github.com/KarypisLab/Gklib
PACKAGEVERSION = git
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "================ Configure"
	@echo "note: custom configure rule"
	@echo "    NOX86 nonzero for non-x86 archs"
NOX86 =
configure :: 
	@if [ ! -z "${SETX}" ] ; then set -x ; fi \
	 && source ${MAKEINCLUDES}/names.sh \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 && setmpicompilers \
	 && ( reportnames \
	 \
	     && rm -rf $${builddir} \
	     && mkdir -p $${builddir} \
	     && cp Makefile.gklib $${srcdir}/Makefile \
	     && if [ ! -z "${NOX86}" -o "${TACC_SYSTEM}" = "macbookair" ] ; then \
	          export nox86="-D NO_X86" ; fi \
	     && cd $${builddir} \
	     && cmdline="cmake \
	            -D CMAKE_INSTALL_PREFIX=$${installdir} \
	            -D CMAKE_C_COMPILER=$${CC} \
	            -D OPENMP=$${ompflag} \
	            -D CMAKE_C_FLAGS=\"-D_POSIX_C_SOURCE=200809L $${nox86}\" \
	            -D CMAKE_POSITION_INDEPENDENT_CODE=ON \
	          $${srcdir}" \
	     && echo "cmake cmdline: $$cmdline" \
	     && echo "    where nox86=<<$${nox86}>> and ompflag=<<$${ompflag}>>" \
	     && eval $$cmdline \
	    ) 2>&1 | tee $${configurelog}

## NO_X86 for mac
## posix 200809L for getline
## posix 199309L for ???

## 	          CONFIG_FLAGS=" -D CMAKE_INSTALL_PREFIX=$${installdir} -D SHARED=1 -D OPENMP=1 -D CMAKE_C_FLAGS=\"-DHAVE_GETLINE -D_POSIX_C_SOURCE=199309L\" "
build ::
	@if [ ! -z "${SETX}" ] ; then set -x ; fi \
	 && source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 && ( cd $${srcdir} \
	    && make V=1 BUILDDIR=$${builddir} \
	    && make BUILDDIR=$${builddir} install \
	    ) 2>&1 | tee $${installlog}
	 @make varsmodule
include ${MAKEINCLUDES}/Make.vars
include ${MAKEINCLUDES}/Make.public

clean ::
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} "" ${PACKAGEBASENAME} \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && cd $${srcdir} && make distclean

GITREPO = https://github.com/KarypisLab/GKlib.git
## git@github.com:KarypisLab/GKlib.git
include ${MAKEINCLUDES}/Make.git

