################################################################
####
#### Makefile for UDUNITS installation
####
################################################################

PACKAGE = UDUNITS
URL = https://www.unidata.ucar.edu/software/udunits/
PACKAGEVERSION = 2.2.28
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
NOLIB = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
default_install : configure build public

TGZURL = https://downloads.unidata.ucar.edu/udunits/${PACKAGEVERSION}/udunits-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
