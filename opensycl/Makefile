################################################################
####
#### Makefile for Opensycl installation
####
#### build with intel only through oneapi
####
################################################################

PACKAGE = opensycl
PACKAGEVERSION = git
MODULES = boost

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
CMAKEFLAGS = \
    -D BOOST_ROOT=${TACC_BOOST_DIR}/include -D Boost_NO_SYSTEM_PATHS=On \
    -D WITH_CUDA_BACKEND=OFF -DWITH_ROCM_BACKEND=OFF -DWITH_LEVEL_ZERO_BACKEND=OFF \
    -D WITH_ACCELERATED_CPU=OFF

GITREPO = https://github.com/OpenSYCL/OpenSYCL.git
include ${MAKEINCLUDES}/Make.git
