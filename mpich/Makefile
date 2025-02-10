################################################################
####
#### Makefile for Mpich installation
####
################################################################

PACKAGE = MPICH
PACKAGEVERSION = 4.3.0
PACKAGEVERSION_DOWNLOAD = 4.3.0rc4
URL = https://www.mpich.org
MODE = seq
FAMILY = mpi
NOMODULE = 1

include ${MAKEINCLUDES}/Make.info
info ::
	@echo "NOTE: this does not generate a module."
	@echo "Create the module file by hand to set new module paths"
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
PKGCONFIG = lib/pkgconfig
MODULENAMEALT = MPI
HASBIN = 1
.PHONY: default_install
default_install : configure build module

TGZURL = http://www.mpich.org/static/downloads/${PACKAGEVERSION_DOWNLOAD}/mpich-${PACKAGEVERSION_DOWNLOAD}.tar.gz
include ${MAKEINCLUDES}/Make.download

SUBMODULE = submodule
## include ${MAKEINCLUDES}/Make.git
