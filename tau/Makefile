################################################################
####
#### Makefile for tau
####
################################################################

PACKAGE = tau
URL = https://www.cs.uoregon.edu/research/tau/home.php
PACKAGEVERSION = 2.33
MODE = mpi
MODULES = pdtoolkit

include ${MAKEINCLUDES}/Make.info

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
	     && ./configure -prefix=$${installdir} \
	            -bfd=download -dwarf=download -unwind=download -iowrapper \
	            -mpi -c++=${TACC_CXX} \
	            -pdt=${TACC_PDTOOLKIT_DIR} \
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
## PKGCONFIGLIB = pkgconfig # seems only obscure packages somewhere deep down
HASBIN = 1
BINDIR = x86_64/bin
NOLIB = 1
include ${MAKEINCLUDES}/Make.vars
build :: varsmodule

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build public

TGZURL = http://tau.uoregon.edu/tau.tgz
include ${MAKEINCLUDES}/Make.download
untar ::
	@source ${MAKEINCLUDES}/names.sh \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 && cd $$srcdir \
	 && echo " .. getting ext" \
	 && wget http://tau.uoregon.edu/ext.tgz && tar fxz ext.tgz
