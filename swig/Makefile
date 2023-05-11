################################################################
####
#### Makefile for swig
####
#### https://www.swig.org/
####
################################################################

PACKAGE = swig
PACKAGEVERSION = 4.1.1
MODE = seq
MODULES = pcre2

CONFIGUREFLAGS=--with-pcre2-prefix=${TACC_PCRE2_DIR}
include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
HASBIN = 1
include ${MAKEINCLUDES}/Make.install

TGZURL = http://prdownloads.sourceforge.net/swig/swig-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
