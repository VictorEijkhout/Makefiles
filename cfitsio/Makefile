################################################################
####
#### Makefile for Cfitsio installation
####
#### https://heasarc.gsfc.nasa.gov/fitsio/
####
################################################################

PACKAGE = cfitsio
PACKAGEVERSION = 4.2.0
MODE = seq
NOMODULES = 

include ${MAKEINCLUDES}/Make.info

PKGCONFIG = lib/pkgconfig
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install

TGZURL = https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
