################################################################
####
#### Makefile for LibCeed installation
####
################################################################

PACKAGE = LIBCEED
ABOUT = Code for Efficient Extensible Discretizations
URL = https://github.com/ceed/libceed
PACKAGEVERSION = git
# 1.0.0 is at rc1
CHMODR = 1
MODULES = petsc
MODE = mpi

include ${MAKEINCLUDES}/Make.info

JCOUNT = 12
.PHONY: configure build
PKGCONFIG = lib/pkgconfig
configure ::
	@echo "there is no separate configure rule"
build :: modules
	@\
	source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 \
	 && ( \
	    source ${MAKEINCLUDES}/compilers.sh \
	     && reportcompilers \
	     && requirenonzero srcdir \
	     && pushd $$srcdir \
	     && cmdline="make --no-print-directory \
	            -B OPT=-g -j ${JCOUNT} V=1 \
	            BUILDDIR=$$builddir/ \
	            prefix=$$installdir \
	            install "\
	     && echo "making: $$cmdline" && eval $$cmdline \
	     && popd && make module \
	     ) 2>&1 | tee $$installlog \
	 && echo "See $$installlog"
include ${MAKEINCLUDES}/Make.vars
default_install :
	@make configure build PACKAGEVERSION=git JCOUNT=${JCOUNT}

####
#### Clean
####

include ${MAKEINCLUDES}/Make.clean

####
#### Download stuff
####

include ${MAKEINCLUDES}/Make.git
GITREPO = https://github.com/CEED/libCEED
