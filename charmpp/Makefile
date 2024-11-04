################################################################
####
#### Makefile for Charm++ installation
####
################################################################

PACKAGE = charmpp
URL = https://charmplusplus.org/ \
    https://github.com/UIUC-PPL/charm
PACKAGEVERSION = 8.0.1

MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig

CMAKEFLAGS = \
    -D NETWORK=mpi \
    -D CMAKE_POSITION_INDEPENDENT_CODE=ON
BUILDSTATICLIBS = 1

# include ${MAKEINCLUDES}/Make.cmake
# include ${MAKEINCLUDES}/Make.cbuild
# info ::
# 	@echo "make default_install : configure build"
# .PHONY: default_install
# default_install : configure build

.PHONY: config-help
config-help :
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 && cd $$srcdir \
	 && ./build --help

.PHONY: default_install
# Vista Intel build fails
# https://github.com/charmplusplus/charm/issues/3851
default_install :
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 && ( \
	    cd $$srcdir \
	     && ./build charm++ \
	            $$( if [ "${TACC_SYSTEM}" = "vista" ] ; then echo mpi-linux-arm8 ; else echo mpi-linux-x86_64 ; fi ) \
	            --with-production --enable-tracing \
	            --install-prefix=$${prefixdir} --force \
	    ) 2>&1 | tee $$installlog \
	 && make --no-print-directory varsmodule
include ${MAKEINCLUDES}/Make.vars
HASBIN=1

TGZURL = https://github.com/charmplusplus/charm/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

