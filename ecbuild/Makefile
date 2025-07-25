################################################################
####
#### Makefile for ecbuild installation
####
#### building block for MAPL
####
################################################################

PACKAGE = ecbuild
URL = https://github.com/ecmwf/ecbuild/tree/develop
PACKAGEVERSION = 3.8.2

MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info
HASBIN = 1
NOLIB = 1
build :
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 \
	 && ( \
	    cd $${srcdir} \
	     && mkdir -p bootstrap && cd bootstrap \
	     && buildline="../bin/ecbuild --prefix=$${prefixdir} .. " \
	     && echo "buildline=$$buildline" \
	     && echo "no ctest" >/dev/null \
	     && make install \
	    )  2>&1 | tee $$installlog
EXTRAINSTALLPATHS = CMAKE_MODULE_PATH=share/ecbuild/cmake
include ${MAKEINCLUDES}/Make.vars

.PHONY: default_install
default_install : build module

TGZURL = https://github.com/ecmwf/ecbuild/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean

