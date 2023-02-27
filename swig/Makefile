################################################################
####
#### Makefile for swig
####
#### https://www.swig.org/
####
################################################################

PACKAGE = swig
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 4.1.1
MODE = seq
MODULES = pcre2

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install

TGZURL = http://prdownloads.sourceforge.net/swig/swig-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
