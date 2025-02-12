################################################################
####
#### Makefile for NETHACK installation
####
################################################################

PACKAGE = NETHACK
PACKAGEVERSION = 3.6.7
MODE = seq

include ${MAKEINCLUDES}/Make.info

# PKGCONFIGLIB = pkgconfig
NOLIB = 1
HASBIN = 1
# include ${MAKEINCLUDES}/Make.cmake
# include ${MAKEINCLUDES}/Make.cbuild
default_install : configure build make_module
configure :
	@source ${MAKEINCLUDES}/names.sh \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	    "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 && cd $$srcdir \
	 && sys/unix/setup.sh sys/unix/hints/linux-x11
build :
	@source ${MAKEINCLUDES}/names.sh \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	    "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 && cd $$srcdir \
	 && make --no-print-directory \
	 && mkdir -p $$installdir/bin \
	 && cp src/nethack $$installdir/bin
include ${MAKEINCLUDES}/Make.vars
make_module :
	@source ${MAKEINCLUDES}/names.sh \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	    "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 && make --no-print-directory varsmodule 

TGZURL = https://github.com/NetHack/NetHack/archive/refs/tags/NetHack-${PACKAGEVERSION}_Released.tar.gz
include ${MAKEINCLUDES}/Make.download
