################################################################
####
#### Makefile for Metis installation
####
#### https://github.com/KarypisLab/METIS
#### THIS IS PROBABLY BROKEN AND NO LONGER MAINTAINED
####
#### https://github.com/KarypisLab/METIS/issues/66
####
################################################################

PACKAGE = metisk
PACKAGEVERSION = 5.2.1
#PACKAGEVERSION = git
MODE = mpi
MODULES = gklib

include ${MAKEINCLUDES}/Make.info

configure :: modules
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
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
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && ( cd $${srcdir} \
	    && echo "\
	        BUILDDIR=$${srcdir}/build \
	        "\
	    && make V=1 \
	    && echo "not doing make     BUILDDIR=$${srcdir}/build install" \
	    ) 2>&1 | tee $${installlog} \
	 && make varsmodule \
	 && ( cd $${srcdir} && make distclean )
.PHONY: default_install
info ::
	@echo "default_install : configure build module"
default_install : configure build module public
include ${MAKEINCLUDES}/Make.vars

GITREPO = git@github.com:KarypisLab/Metis.git
include ${MAKEINCLUDES}/Make.git

TGZURL = https://github.com/KarypisLab/METIS/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
