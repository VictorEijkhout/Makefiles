################################################################
####
#### Makefile for rangev3
####
#### https://github.com/ericniebler/range-v3
####
################################################################

PACKAGE = rangev3
PACKAGEVERSION = 0.12.0
MODE = seq

include ${MAKEINCLUDES}/Make.info

# intel compiler considers missing proto as fatal
CMAKEFLAGS = --compile-no-warning-as-error
# -DCMAKE_COMPILE_WARNING_AS_ERROR=NO
# -DCMAKE_CXX_FLAGS=-Wno-error
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild


TGZURL = https://github.com/ericniebler/range-v3/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
