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
	     && make BUILDDIR=$${builddir} prefix=$${installdir} \
	            config shared=1 openmp=ON \
	            cc-flags="-D_POSIX_C_SOURCE=199309L" \
	    ) 2>&1 | tee $${configurelog}
build ::
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} "" ${PACKAGEBASENAME} \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && ( cd $${srcdir} \
	    && make BUILDDIR=$${builddir} \
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

