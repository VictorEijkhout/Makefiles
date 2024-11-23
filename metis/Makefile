################################################################
####
#### Makefile for Metis as patched by scivision installation
####
#### https://github.com/scivision/METIS
#### this the 5.1.0.3 version only
#### there are big problems with the Karypis 5.2 version
####
####   --with-64-bit-indices option requires a metis build with IDXTYPEWIDTH=64
#### however
#### https://github.com/scivision/METIS/issues/7
####
################################################################

PACKAGE = metis
PACKAGEVERSION = 5.1.0.3
URL = https://github.com/scivision/METIS
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info
HASBIN = 1
include ${MAKEINCLUDES}/Make.cmake
WIDTH=32
CMAKEFLAGS = -DCMAKE_C_FLAGS='-DIDXTYPEWIDTH=${WIDTH} -DREALTYPEWIDTH=${WIDTH}'
include ${MAKEINCLUDES}/Make.cbuild

.PHONY: default_install i32 i64 i64-configure i64-build
info ::
	@echo "make default_install : i32; also i64 target"
default_install : i32
i32 : 
	@make --no-print-directory \
	    configure build public \
	    INSTALLVARIANT=i32 
i64 : i64-configure i64-build
# stupid post-processing:
# https://github.com/scivision/METIS/issues/7
i64-configure :
	@make --no-print-directory \
	    configure \
	    WIDTH=64 \
	    INSTALLVARIANT=i64 MODULEVERSIONEXTRA=i64
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}"  "${PACKAGEVERSION}"  "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}"  "${INSTALLROOT}"  "${INSTALLEXT}"  "i64" \
	    "${HOMEDIR}"  "${BUILDDIRROOT}"  "${MODE}" \
	    "${PREFIXOPTION}"  "${PREFIXEXTRA}" \
	 && sed -i $${builddir}/src/include/metis.h \
	    -e '/define/s/WIDTH 32/WIDTH 64/'
i64-build :
	@make --no-print-directory \
	    build public \
	    WIDTH=64 \
	    INSTALLVARIANT=i64 MODULEVERSIONEXTRA=i64

TGZURL = https://github.com/scivision/METIS/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

