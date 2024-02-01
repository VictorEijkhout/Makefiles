################################################################
####
#### Makefile for Boost installation
####
################################################################

PACKAGE = boost
PACKAGEVERSION = 1.83.0
PACKAGEVERSION_DOWNLOAD = 1_83_0
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
	     && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "" "" \
	 && ( \
	    rm -rf $$builddir && mkdir -p $$builddir \
	     && mkdir -p $$installdir && rm -rf $$installdir/* \
	     && cd $$srcdir \
	     && ./bootstrap.sh --prefix=$$installdir \
	        --with-libraries=atomic,chrono,container,context,contract,coroutine,date_time,exception,fiber,filesystem,graph,graph_parallel,headers,iostreams,json,locale,log,math${MPI},nowide,program_options,python,random,regex,serialization,stacktrace,system,test,thread,timer,type_erasure,url,wave \
	     && ./b2 toolset=gcc \
	          variant=release link=shared debug-symbols=off inlining=full \
	          address-model=64  stage --with-fiber --with-context --with-test --with-atomic --with-filesystem --with-headers install \
	          --prefix=$${installdir} --build-dir=$${builddir} \
	     && echo "not /b2 --with-fiber --with-context --with-filesystem --with-atomic --with-test" \
	     && echo "not ./b2 install --prefix=$$installdir --build-dir=$$builddir toolset=gcc" \
	     && echo "inoperative option: link=shared" >/dev/null \
	    ) 2>&1 | tee $$installlog \
	 && make --no-print-directory varsmodule \
	      CMAKEPREFIXLIB=cmake
.PHONY: seq par default_install
default_install : seq
seq : configure build MPI=
par : configure build MPI=,mpi

include ${MAKEINCLUDES}/Make.vars
## NOLIB = 1
include ${MAKEINCLUDES}/Make.public

TGZURL = https://boostorg.jfrog.io/artifactory/main/release/${PACKAGEVERSION}/source/boost_${PACKAGEVERSION_DOWNLOAD}.tar.gz
include ${MAKEINCLUDES}/Make.download
