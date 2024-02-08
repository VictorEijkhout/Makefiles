################################################################
####
#### Makefile for Adaptivecpp installation
####
################################################################

PACKAGE = adaptivecpp
URL = https://github.com/AdaptiveCpp/AdaptiveCpp
PACKAGEVERSION = 23.10.0
MODULES = boost clang
MODE = seq

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
HASBIN = 1
include ${MAKEINCLUDES}/Make.cbuild
HASCUDA = OFF
CMAKEFLAGS = \
    $$( if [ ${PACKAGEVERSION} != "git" ] ; then echo '-D ACPP_VERSION_SUFFIX="-"' ; fi ) \
    -D BOOST_ROOT=${TACC_BOOST_DIR} -D Boost_NO_SYSTEM_PATHS=On -DBoost_DEBUG=ON \
    -D WITH_CUDA_BACKEND=${HASCUDA} -DWITH_ROCM_BACKEND=OFF -DWITH_LEVEL_ZERO_BACKEND=OFF \
    \
    -D BUILD_CLANG_PLUGIN=OFF \
    -D CMAKE_CXX_COMPILER=${TACC_CXX} \
    -D WITH_ACCELERATED_CPU=OFF \
    \
    -D FIXES_FOR_MISSING_LIB_CLANG=ON \
    -D WITH_SSCP_COMPILER=OFF \
    -D WITH_OPENCL_BACKEND=OFF \
    -D WITH_LEVEL_ZERO_BACKEND=OFF

.PHONY: cpu gpu
info ::
	@echo "make cpu / gpu (load cuda or not)"
cpu : 
	@make modules
	@make --no-print-directory HASCUDA=OFF NONMODULES=cuda \
	    configure build
gpu : 
	@if [ "${TACC_SYSTEM}" = "frontera" ] ; then \
	    echo && echo "Last time I looked Frtonera did not have LLVM" && echo \
	     && exit 1 ; fi
	@make --no-print-directory HASCUDA=ON \
	    configure build
.PHONY: default_install
info ::
	@echo "make default_install : cpu version"
default_install : cpu
	@echo && echo "==== cpu install done, not gpu ====" && echo

TGZURL = https://github.com/AdaptiveCpp/AdaptiveCpp/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/Adaptivecpp/Adaptivecpp.git
include ${MAKEINCLUDES}/Make.git
