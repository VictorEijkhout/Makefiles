################################################################
####
#### Makefile for MUMPS installation
####
################################################################

PACKAGE = MUMPS
PACKAGEVERSION = 5.8.0
URL = https://www.mumps.org/
MODE = mpi

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

TGZURL = https://mumps-solver.org/MUMPS_${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
