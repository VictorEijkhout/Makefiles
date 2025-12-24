################################################################
####
#### Makefile for MPUNITS installation
####
################################################################

PACKAGE = MPUNITS
ABOUT = The quantities and units library for C++
PACKAGEVERSION = 2.5.0
URL = https://github.com/mpusz/mp-units
MODE = seq
MODULES = fmtlib gsllite
# only needed for the full library install: catch2

include ${MAKEINCLUDES}/Make.info

CMAKEPREFIXPATHSET = 1
NOLIB = 1

# there is a top level cmake setup but that bombs in all sorts of ways.
CMAKESUBDIR = src

CMAKEFLAGS = \
    -D MP_UNITS_API_STD_FORMAT=OFF \
    -D MP_UNITS_LIBCXX=ON

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build public

DOWNLOADURL = https://github.com/mpusz/mp-units/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/mpusz/mp-units.git
include ${MAKEINCLUDES}/Make.git

    # if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "2024")
    #     message(STATUS "Intel 1API 2024+ with libc++ detected, overriding `std::format` support")
    #     set(${projectPrefix}LIB_FORMAT_SUPPORTED ON)
    # endif()

# gcc compile unit test fail:
# [ 63%] ESC[32mESC[1mLinking CXX executable unit_tests_runtimeESC[0m
# cd /work2/00434/eijkhout/mpunits/build-2.4.0-stampede3-gcc13.2/test/runtime && /home1/apps/cmake/3.31.5/bin/cmake -E cmake_link_script CMakeFiles/unit_tests_runtime.dir/link.txt --verbose=1
# /opt/apps/gcc/13.2.0/bin/ld: /work2/00434/eijkhout/catch2/installation-catch2-3.5.4-stampede3-gcc13.2/lib64/libCatch2.so.3.5.4: undefined reference to `__cxa_call_terminate@CXXABI_1.3.15'
# collect2: error: ld returned 1 exit status
