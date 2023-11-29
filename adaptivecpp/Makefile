################################################################
####
#### Makefile for Adaptivecpp installation
####
#### build with intel only through oneapi
####
################################################################

PACKAGE = adaptivecpp
URL = https://github.com/AdaptiveCpp/AdaptiveCpp
PACKAGEVERSION = 23.10.0
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

.PHONY: cpu gpu default_install
cpu : 
	@make --no-print-directory HASCUDA=OFF \
	    configure build
gpu : 
	@make --no-print-directory HASCUDA=ON  \
	    configure build
.PHONY: default_install
info ::
	@echo "make default_install : cpu version"

TGZURL = https://github.com/AdaptiveCpp/AdaptiveCpp/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/Adaptivecpp/Adaptivecpp.git
include ${MAKEINCLUDES}/Make.git
