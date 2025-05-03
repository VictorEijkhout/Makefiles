################################################################
####
#### Makefile for swig
####
################################################################

PACKAGE = swig
URL = https://www.swig.org/ https://github.com/swig/swig
PACKAGEVERSION = 4.3.1
MODE = seq
MODULES = pcre2 bison

include ${MAKEINCLUDES}/Make.info

# cmake install gives error:
# https://github.com/swig/swig/issues/3171

# # Stupid: Swig has its own FindPCRE2 rather than relying on cmake:
# # swig-4.1.1/Tools/cmake/FindPCRE2.cmake
# CMAKEFLAGS = \
#     -D PCRE2_INCLUDE_DIR:PATH=${TACC_PCRE2_INC} \
#     -D PCRE2_LIBRARY:PATH=${TACC_PCRE2_LIB}/libpcre2-8.so
# include ${MAKEINCLUDES}/Make.cmake
# include ${MAKEINCLUDES}/Make.cbuild

include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
.PHONY: default_install
default_install : configure build
HASBIN = 1
NOLIB = 1

TGZURL = http://prdownloads.sourceforge.net/swig/swig-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
