################################################################
####
#### Makefile for OSUbenchmark installation
####
#### https://github.com/icecube/OSUbenchmark
####
################################################################

PACKAGE = OSUBENCHMARK
PACKAGEVERSION = 6.0
MODULES =
MODE = mpi

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
NOLIB = 1
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

TGZURL = http://mvapich.cse.ohio-state.edu/download/mvapich/osu-micro-benchmarks-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
