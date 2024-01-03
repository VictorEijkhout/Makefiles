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
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install : configure build

TGZURL = http://prdownloads.sourceforge.net/swig/swig-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
