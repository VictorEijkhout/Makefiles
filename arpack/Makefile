################################################################
####
#### Makefile for Arpack installation
#### 
################################################################

PACKAGE = ARPACK
URL = https://github.com/opencollab/arpack-ng
PACKAGEVERSION = 3.9.1
MODE = seq
MODULES = eigen

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
    -D EXAMPLES=ON -D ICB=ON -D EIGEN=ON

include ${MAKEINCLUDES}/Make.cmake
# HASBIN=1
# PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make par : parallel installation"
	@echo "make seq : sequential installation"
ARPACKPARALLEL = OFF
.PHONY: seq par
seq :
	@make configure build JCOUNT=${JCOUNT} \
	    PACKAGEVERSION=${PACKAGEVERSION} \
	    ARPACKPARALLEL=OFF MODE=seq MODULENAME=arpack
par :
	@make configure build JCOUNT=${JCOUNT} \
	    PACKAGEVERSION=${PACKAGEVERSION} \
	    ARPACKPARALLEL=ON  MODE=mpi MODULENAME=parpack MODULENAMEALT=arpack
info ::
	@echo "make default_install : sequential and parallel versions"
default_install :: 
	@make --no-print-directory seq par JCOUNT=${JCOUNT} PACKAGEVERSION=${PACKAGEVERSION} \
	    seq par

##
## download
##
TGZURL = https://github.com/opencollab/arpack-ng/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

##
## git clone
##
GITREPO = https://github.com/opencollab/arpack-ng.git
include ${MAKEINCLUDES}/Make.git

##
## Clean 
##
include ${MAKEINCLUDES}/Make.clean
