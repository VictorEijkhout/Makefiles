################################################################
####
#### Makefile for Alps installation
####
################################################################

PACKAGE = ALPS
URL = http://alps.space/ \
    https://github.com/danielver02/ALPS
PACKAGEVERSION = 1.0.1
MODE = mpi
MODULES = mkl

include ${MAKEINCLUDES}/Make.info
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
HASBIN=1

.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build

TGZURL = https://github.com/danielver02/ALPS/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

