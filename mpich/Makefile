################################################################
####
#### Makefile for Mpich installation
####
#### https://www.mpich.org/2022/01/12/mpich-4-0rc3-released/
####
################################################################

PACKAGE = MPICH
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 4.1
MODE = seq
FAMILY = mpi

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
PKGCONFIGSET = lib/pkgconfig

TGZURL = http://www.mpich.org/static/downloads/${PACKAGEVERSION}/mpich-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

SUBMODULE = submodule
## include ${MAKEINCLUDES}/Make.git
