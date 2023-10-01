################################################################
####
#### Makefile for Opensycl installation
####
#### build with intel only through oneapi
####
################################################################

PACKAGE = opensycl
URL = https://github.com/OpenSYCL/OpenSYCL
PACKAGEVERSION = 0.9.4
MODULES = boost
MODE = seq

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
HASBIN = 1
include ${MAKEINCLUDES}/Make.cbuild
HASCUDA = OFF
CMAKEFLAGS = \
    -D BOOST_ROOT=${TACC_BOOST_DIR} -D Boost_NO_SYSTEM_PATHS=On \
    -D WITH_CUDA_BACKEND=${HASCUDA} -DWITH_ROCM_BACKEND=OFF -DWITH_LEVEL_ZERO_BACKEND=OFF \
    -D WITH_ACCELERATED_CPU=OFF

.PHONY: cpu gpu
cpu : 
	@make --no-print-directory HASCUDA=OFF \
	    configure build
gpu : 
	@make --no-print-directory testmodule MODULE=cuda
	@make --no-print-directory HASCUDA=ON  \
	    configure build

TGZURL = https://github.com/OpenSYCL/OpenSYCL/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/OpenSYCL/OpenSYCL.git
include ${MAKEINCLUDES}/Make.git
