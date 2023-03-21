################################################################
####
#### Makefile for Boost installation
####
################################################################

PACKAGE = boost
PACKAGEVERSION = 1_81_0
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "================ bootstrap install"
	@echo "make build"
.PHONY: configure build
configure ::
	@echo "there is no configure rule"
build ::
	@source ${MAKEINCLUDES}/names.sh \
	     && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	     && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} ${INSTALLEXT} \
	 && ( \
	    rm -rf $$builddir && mkdir -p $$builddir \
	     && rm -rf $$installdir && mkdir -p $$installdir \
	     && cd $$srcdir \
	     && ./bootstrap.sh --prefix=$$installdir \
	        --with-libraries=atomic,chrono,container,context,contract,coroutine,date_time,exception,fiber,filesystem,graph,graph_parallel,headers,iostreams,json,locale,log,math,mpi,nowide,program_options,python,random,regex,serialization,stacktrace,system,test,thread,timer,type_erasure,url,wave \
	     && ./b2 install --prefix=$$installdir --build-dir=$$builddir --toolset=intel \
	    ) 2>&1 | tee $@.log \
	 && make --no-print-directory varsmodule
include ${MAKEINCLUDES}/Make.vars
include ${MAKEINCLUDES}/Make.public

TGZURL = https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
