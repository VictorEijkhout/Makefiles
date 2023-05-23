################################################################
####
#### Makefile for Netcdf installation
####
################################################################

PACKAGE = SUITESPARSE
PACKAGEVERSION = 7.0.1
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

## include ${MAKEINCLUDES}/Make.cmake
## include ${MAKEINCLUDES}/Make.cbuild

.PHONY: configure build
configure :
	@echo "There is no separate configure rule"
build :
	@source ${MAKEINCLUDES}/names.sh \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 && setmpicompilers \
	 && reportcompilers && echo \
	 && ( \
	    cd $$srcdir \
	     && CMAKE_OPTIONS="\
	        -D CMAKE_C_COMPILER=$${cc} -D CMAKE_CXX_COMPILER=$${cxx} \
	        -D CMAKE_INSTALL_PREFIX=$$installdir \
	                      " \
	        make \
	     && make install \
	    ) 2>&1 | tee $$installlog \
	 && echo && echo "Installation finished in: $$installdir" && echo


TGZURL = https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

