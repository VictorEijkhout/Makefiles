################################################################
####
#### Makefile for tau
####
################################################################

PACKAGE = tau
URL = https://www.cs.uoregon.edu/research/tau/home.php
PACKAGEVERSION = 2.33
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

HASBIN = 1
NOLIB = 1
# include ${MAKEINCLUDES}/Make.cmake
# include ${MAKEINCLUDES}/Make.cbuild
.PHONY: configure build
info ::
	@echo "make configure build"
configure ::
	@source ${MAKEINCLUDES}/names.sh \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 && ( \
	    cd $$srcdir \
	     && ./configure \
	            -bfd=download -dwarf=download -unwind=download -iowrapper \
	            -pdt=pdtoolkit-3.25.1 \
	    ) \
	    | tee $${configurelog}
build :: modules
	@source ${MAKEINCLUDES}/names.sh \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 && setmpicompilers && reportcompilers \
	 \
	 && ( \
	    cd $$srcdir \
	     && make install \
	    ) 2>&1 | tee $$installlog \
	 && echo "See $$installlog"

TGZURL = http://tau.uoregon.edu/tau.tgz
include ${MAKEINCLUDES}/Make.download
untar ::
	@source ${MAKEINCLUDES}/names.sh \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 && cd $$srcdir \
	 && echo " .. getting ext" \
	 && wget http://tau.uoregon.edu/ext.tgz && tar fxz ext.tgz \
	 && echo " .. getting pd5" \
	 && wget http://tau.uoregon.edu/pdt.tgz && tar fxz pdt.tgz
