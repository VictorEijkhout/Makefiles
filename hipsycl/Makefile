################################################################
####
#### Makefile for Hipsycl installation
####
################################################################

PACKAGE = hipsycl
PACKAGEVERSION = git
MODULES = gcc boost

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
CMAKEFLAGS = -D BOOST_ROOT=${TACC_BOOST_DIR} -D Boost_NO_SYSTEM_PATHS=On \
    -D WITH_CUDA_BACKEND=OFF -DWITH_ROCM_BACKEND=OFF -DWITH_LEVEL_ZERO_BACKEND=OFF \
    -D WITH_ACCELERATED_CPU=OFF

GITREPO = https://github.com/illuhad/hipSYCL.git
include ${MAKEINCLUDES}/Make.git
