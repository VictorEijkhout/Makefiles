################################################################
####
#### Makefile for Metis as patched by scivision installation
####
#### https://github.com/scivision/METIS
#### this the 5.1.0.3 version only
#### there are big problems with the Karypis 5.2 version
####
####   --with-64-bit-indices option requires a metis build with IDXTYPEWIDTH=64
#### however
#### https://github.com/scivision/METIS/issues/7
####
#### https://github.com/scivision/METIS/issues/10
#### https://github.com/scivision/METIS/issues/8
####
################################################################

PACKAGE = metis
ABOUT = CMake enhanced Karypis METIS library
URL = https://github.com/scivision/METIS/

# PACKAGEVERSION = 5.2.1.2
# still problems with 5.2.1.2:
## https://github.com/scivision/METIS/issues/8

PACKAGEVERSION=git
GITDATE=20250821

MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info
include ${MAKEINCLUDES}/Make.cmake
WIDTH=32
# 1. posix source:
#    https://github.com/KarypisLab/GKlib/issues/24
# 2. testing off to prevent sqrt@@GLIBC error
# 3. copy-dt-needed-entries because 
#    build-5.2.1.1-clx-gcc13.2-impi21.9.0-i32/_deps/gklib-build/libGKlib.so
#    see https://stackoverflow.com/a/55086637/2044454
ifeq "${TACC_FAMILY_COMPILER}" "nvidia"
  EXTRA_C_FLAGS = -D_POSIX_C_SOURCE=200809L
else
  EXTRA_C_FLAGS = -D_POSIX_C_SOURCE=200809L -Wl,--copy-dt-needed-entries
endif
CMAKEFLAGS = \
    -D UNIX=ON \
    -D BUILD_TESTING=OFF \
    -D IDXTYPEWIDTH=${WIDTH} -D REALTYPEWIDTH=${WIDTH} \
    -D CMAKE_C_FLAGS='${EXTRA_C_FLAGS}'
include ${MAKEINCLUDES}/Make.cbuild

.PHONY: default_install i32 i64 i64-configure i64-build
info ::
	@echo "make default_install : i32; also i64 target"
HASBIN = 1
default_install : i32
i32 : 
	@make --no-print-directory \
	    configure build public \
	    INSTALLVARIANT=i32 
i64 : i64-configure i64-build
i64-configure :
	@make --no-print-directory \
	    configure \
	    WIDTH=64 \
	    INSTALLVARIANT=i64 MODULEVERSIONEXTRA=i64
i64-build :
	@make --no-print-directory \
	    build public \
	    WIDTH=64 \
	    INSTALLVARIANT=i64 MODULEVERSIONEXTRA=i64
# stupid post-processing is solved in the 5.2.1.2 release:
# https://github.com/scivision/METIS/issues/7
i64-fix :
	/bin/true \
	 && sed -i $${builddir}/src/include/metis.h \
	    -e '/define/s/WIDTH 32/WIDTH 64/'

TGZURL = https://github.com/scivision/METIS/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
# BOOTSTRAP = sed -i \
#     -e '/GIT_TAG/s/6897bfe/078042614ff5192ab4c9bde0575bfef0d4572687/' \
#     -e '/GIT_TAG/s/44630ca/078042614ff5192ab4c9bde0575bfef0d4572687/' \
#     CMakeLists.txt
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/scivision/METIS.git
include ${MAKEINCLUDES}/Make.git

ifdef TICKETS
intel24 compilation S3
/opt/apps/gcc/13.2.0/bin/ld: CMakeFiles/gpmetis.dir/gpmetis.c.o: undefined reference to symbol 'gk_GetCurMemoryUsed'
/opt/apps/gcc/13.2.0/bin/ld: /work2/00434/eijkhout/metis/build-5.2.1.1-stampede3-intel24.2-impi21.13-i32/_deps/gklib-build/libGKlib.so: error adding symbols: DSO missing from command line

Attempt to use our own gklib in CMakeLists.txt
This does not work. Really needs the included gklib
## FetchContent_Declare(GKlib
## GIT_REPOSITORY https://github.com/scivision/GKlib.git
## #GIT_TAG 44630ca
## GIT_SHALLOW true
## )
## FetchContent_MakeAvailable(GKlib)
find_package( GKlib REQUIRED ) # VLE

endif
