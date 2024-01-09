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
# Stupid: Swig has its own FindPCRE2 rather than relying on cmake:
# swig-4.1.1/Tools/cmake/FindPCRE2.cmake
CMAKEFLAGS = \
    -D PCRE2_INCLUDE_DIR:PATH=${TACC_PCRE2_INC} \
    -D PCRE2_LIBRARY:PATH=${TACC_PCRE2_LIB}/libpcre2-8.so
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install : configure build

TGZURL = http://prdownloads.sourceforge.net/swig/swig-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
