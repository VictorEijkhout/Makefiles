################################################################
####
#### Makefile for swig
####
################################################################

PACKAGE = swig
URL = https://www.swig.org/
PACKAGEVERSION = 4.1.1
MODE = seq
MODULES = pcre2

include ${MAKEINCLUDES}/Make.info

HASBIN = 1
NOLIB = 1
# CONFIGUREFLAGS=--with-pcre2-prefix=${TACC_PCRE2_DIR}
# include ${MAKEINCLUDES}/Make.configure
# include ${MAKEINCLUDES}/Make.install
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = http://prdownloads.sourceforge.net/swig/swig-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
