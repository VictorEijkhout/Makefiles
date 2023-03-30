################################################################
####
#### Makefile for GKlib installation
####
#### https://github.com/KarypisLab/Gklib
####
################################################################

PACKAGE = gklib
PACKAGEVERSION = git
MODE = mpi
MODULES = gklib metis

include ${MAKEINCLUDES}/Make.info

configure :: 
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} "" ${PACKAGEBASENAME} \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 \
	 && rm -rf $${builddir} $${installdir} \
	 && cp Makefile.gklib $${srcdir}/Makefile \
	 && ( cd $${srcdir} \
	     && make BUILDDIR=$${builddir} \
	          CONFIG_FLAGS=" -D CMAKE_INSTALL_PREFIX=$${installdir} -D CMAKE_C_FLAGS=\"-DHAVE_GETLINE -D_POSIX_C_SOURCE=199309L\" " \
	            config \
	    ) 2>&1 | tee $${configurelog}
## 	          CONFIG_FLAGS=" -D CMAKE_INSTALL_PREFIX=$${installdir} -D SHARED=1 -D OPENMP=1 -D CMAKE_C_FLAGS=\"-DHAVE_GETLINE -D_POSIX_C_SOURCE=199309L\" "
build ::
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} "" ${PACKAGEBASENAME} \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && ( cd $${srcdir} \
	    && make V=1 BUILDDIR=$${builddir} \
	    && make BUILDDIR=$${builddir} install \
	    ) 2>&1 | tee $${installlog}
	 @make varsmodule
include ${MAKEINCLUDES}/Make.vars

clean ::
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} "" ${PACKAGEBASENAME} \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && cd $${srcdir} && make distclean

GITREPO = git@github.com:KarypisLab/GKlib.git
include ${MAKEINCLUDES}/Make.git

