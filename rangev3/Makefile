################################################################
####
#### Makefile for rangev3
####
################################################################

URL = https://github.com/ericniebler/range-v3
PACKAGE = rangev3
PACKAGEVERSION = 0.12.0
MODE = seq

include ${MAKEINCLUDES}/Make.info

# intel compiler considers missing proto as fatal
CMAKEFLAGS = \
    -DCMAKE_COMPILE_WARNING_AS_ERROR=OFF
##    -DCMAKE_CXX_COMPILER_FLAGS=-Wno-error
nocmakeflag = \
    --compile-no-warning-as-error \
    -DCMAKE_CXX_FLAGS=-Wno-error
include ${MAKEINCLUDES}/Make.cmake
fix :
	source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${\PACKAGEBASENAME}" "${VARIANT}" \
	 && LANG=C find $$builddir \
	      -type f \
	      -exec sed -i -e 's/Werror/Wno-error/g' {} \; 

include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/ericniebler/range-v3/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/ericniebler/range-v3.git
include ${MAKEINCLUDES}/Make.git
