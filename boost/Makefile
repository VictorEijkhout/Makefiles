################################################################
####
#### Makefile for Boost installation
####
################################################################

PACKAGE = boost
PACKAGEVERSION = 1.84.0
PACKAGEVERSION_DOWNLOAD = 1_84_0
URL = https://www.boost.org/ \
    https://www.boost.org/doc/libs/${PACKAGEVERSION_DOWNLOAD}/tools/build/doc/html/index.html
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "================ bootstrap install"
	@echo "make build"
.PHONY: configure build
configure :: modules
	@echo "there is no configure rule"
build ::
	@source ${MAKEINCLUDES}/names.sh \
	     && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	     && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	            "${PACKAGEBASENAME}" "" "" \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 && ( \
	    rm -rf $$builddir && mkdir -p $$builddir \
	     && mkdir -p $$installdir && rm -rf $$installdir/* \
	     && cd $$srcdir \
	     && ./bootstrap.sh --prefix=$$installdir \
	     && ./b2 toolset=${TACC_FAMILY_COMPILER} \
	          variant=release link=shared debug-symbols=off inlining=full \
	          address-model=64  stage \
	            --with-atomic --with-chrono --with-container --with-context \
	            --with-contract --with-coroutine --with-date_time \
	            --with-exception --with-fiber --with-filesystem --with-graph --with-graph_parallel \
	            --with-headers --with-iostreams --with-json --with-locale --with-log \
	            --with-math --with-nowide --with-program_options \
	            --with-random --with-regex --with-serialization --with-stacktrace --with-system \
	            --with-test --with-thread --with-timer \
	            --with-type_erasure --with-url --with-wave \
	            install \
	            --prefix=$${installdir} --build-dir=$${builddir} \
	     && echo "inoperative option: link=shared" >/dev/null \
	    ) 2>&1 | tee $$installlog \
	 && make --no-print-directory varsmodule \
	      CMAKEPREFIXPATHSET=1

.PHONY: seq par default_install
default_install : seq
seq :
	@make --no-print-directory configure build \
	    PACKAGEVERSION=${PACKAGEVERSION} JCOUNT=${JCOUNT} \
	    MODE=seq MODULENAME=boost MPI= 
par : 
	@make --no-print-directory configure build \
	    PACKAGEVERSION=${PACKAGEVERSION} JCOUNT=${JCOUNT} \
	    MODE=mpi MODULENAME=boostmpi MPI=,mpi

include ${MAKEINCLUDES}/Make.vars
## NOLIB = 1
include ${MAKEINCLUDES}/Make.public

TGZURL = https://boostorg.jfrog.io/artifactory/main/release/${PACKAGEVERSION}/source/boost_${PACKAGEVERSION_DOWNLOAD}.tar.gz
include ${MAKEINCLUDES}/Make.download
