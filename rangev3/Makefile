################################################################
####
#### Makefile for rangev3
####
#### Intel problems:
#### https://github.com/ericniebler/range-v3/issues/1845
####
#### Vista gcc aa64 problem:
#### https://github.com/ericniebler/range-v3/issues/1846
####
################################################################

URL = https://github.com/ericniebler/range-v3
PACKAGE = rangev3
PACKAGEVERSION = 0.12.0
NONMODULES = clang
MODE = seq

include ${MAKEINCLUDES}/Make.info

# intel compiler considers missing proto as fatal
CMAKEFLAGS = \
    -D RANGES_ENABLE_WERROR=OFF \
    -DCMAKE_COMPILE_WARNING_AS_ERROR=OFF
include ${MAKEINCLUDES}/Make.cmake
fix :
	source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${\PACKAGEBASENAME}" "${VARIANT}" \
	 && LANG=C find $$builddir \
	      -type f \
	      -exec sed -i -e 's/-Werror//gp' {} \; 

CMAKEPREFIXLIB = cmake/range-v3
include ${MAKEINCLUDES}/Make.cbuild

.PHONY: default_install
default_install : configure build

DOWNLOADURL = https://github.com/ericniebler/range-v3/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/ericniebler/range-v3.git
include ${MAKEINCLUDES}/Make.git
