################################################################
####
#### Makefile for Cfitsio installation
####
#### https://heasarc.gsfc.nasa.gov/fitsio/
####
################################################################

PACKAGE = cfitsio
ABOUT = C and Fortran library for reading and writing FITS files
URL = https://github.com/HEASARC/cfitsio
PACKAGEVERSION = 4.6.2
MODE = seq
NOMODULES = 

include ${MAKEINCLUDES}/Make.info

PKGCONFIG = lib/pkgconfig
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build public

TGZURL = https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
