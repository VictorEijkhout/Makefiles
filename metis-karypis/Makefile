################################################################
####
#### Makefile for Metis installation
####
#### https://github.com/KarypisLab/METIS
#### THIS IS PROBABLY BROKEN AND NO LONGER MAINTAINED
####
################################################################

PACKAGE = metis
PACKAGEVERSION = 5.2.1
#PACKAGEVERSION = git
MODE = mpi
MODULES = gklib

include ${MAKEINCLUDES}/Make.info

configure :: modules
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} "" ${PACKAGEBASENAME} \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 && setmpicompilers \
	 \
	 && requirenonzero builddir \
	 && rm -rf $${builddir} && mkdir -p $${builddir} \
	 && echo "cp Makefile.metis   $${srcdir}/Makefile "\
	 && echo "cp CMakeLists.metis $${srcdir}/CMakeLists.txt "\
	 && ( cd $${srcdir} \
	     && requirenonzero installdir \
	     && reportcompilers \
	     && echo "\
	            BUILDDIR=$${srcdir}/build \
	        "\
	     && make \
	            CONFIG_FLAGS="\
	                -DCMAKE_VERBOSE_MAKEFILE=1 \
	                -D CMAKE_INSTALL_PREFIX=$${installdir} \
	                -D CMAKE_C_COMPILER=$${CC} \
	                -D OPENMP=$${ompflag} -D SHARED=1 \
	                -D CMAKE_POSITION_INDEPENDENT_CODE=ON \
	                " \
	            config \
	    ) 2>&1 | tee $${configurelog}
## not quite working BUILDDIR=$${builddir}, see distclean below
build ::
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} "" ${PACKAGEBASENAME} \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && ( cd $${srcdir} \
	    && echo "\
	        BUILDDIR=$${srcdir}/build \
	        "\
	    && make V=1 \
	    && echo "not doing make     BUILDDIR=$${srcdir}/build install" \
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

TGZURL = https://github.com/KarypisLab/METIS/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
