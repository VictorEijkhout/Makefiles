################################################################
####
#### Makefile for OSUbenchmark installation
####
#### https://github.com/icecube/OSUbenchmark
####
################################################################

PACKAGE = OSUBENCHMARK
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 6.0
MODULES =
MODE = mpi

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install

TGZURL = http://mvapich.cse.ohio-state.edu/download/mvapich/osu-micro-benchmarks-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
