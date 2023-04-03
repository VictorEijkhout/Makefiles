################################################################
####
#### Makefile for Metis installation
####
#### https://github.com/KarypisLab/METIS
####
################################################################

PACKAGE = metis
PACKAGEVERSION = git
MODE = mpi
MODULES = gklib

include ${MAKEINCLUDES}/Make.info

configure :: 
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} "" ${PACKAGEBASENAME} \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 && setmpicompilers \
	 && reportcompilers \
	 \
	 && requirenonzero builddir \
	 && rm -rf $${builddir} && mkdir -p $${builddir} \
	 && cp Makefile.metis   $${srcdir}/Makefile \
	 && cp CMakeLists.metis $${srcdir}/CMakeLists.txt \
	 && ( cd $${srcdir} \
	     && requirenonzero installdir \
	     && make config \
	            BUILDDIR=$${srcdir}/build \
	            prefix=$${installdir} \
	            cc=$${CC} shared=1 openmp=ON \
	            gklib_path=${TACC_GKLIB_DIR} \
	    ) 2>&1 | tee $${configurelog}
## not quite working BUILDDIR=$${builddir}, see distclean below
build ::
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} "" ${PACKAGEBASENAME} \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && ( cd $${srcdir} \
	    && make V=1 BUILDDIR=$${srcdir}/build \
	    && make     BUILDDIR=$${srcdir}/build install \
	    ) 2>&1 | tee $${installlog} \
	 && make varsmodule \
	 && ( cd $${srcdir} && make distclean )
include ${MAKEINCLUDES}/Make.vars

clean ::
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} "" ${PACKAGEBASENAME} \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && cd $${srcdir} && make distclean

GITREPO = git@github.com:KarypisLab/Metis.git
include ${MAKEINCLUDES}/Make.git

