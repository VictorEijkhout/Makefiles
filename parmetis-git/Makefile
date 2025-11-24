################################################################
####
#### Makefile for Parmetis installation
####
#### http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download
####
#### this compiles on my laptop but not stampede2
####
################################################################

PACKAGE = parmetis
PACKAGEVERSION = 4.0.3
MODE = mpi

include ${MAKEINCLUDES}/Make.info

configure :: 
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} "" ${PACKAGEBASENAME} \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 \
	 && rm -rf $${builddir} $${installdir} \
	 && cp Makefile.parmetis $${srcdir} \
	 && ( cd $${srcdir} \
	     && make BUILDDIR=$${builddir} prefix=$${installdir} \
	            config shared=1 openmp=ON \
	            cc-flags=-fPIC cxx-flags=-fPIC \
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

DOWNLOADURL = http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

