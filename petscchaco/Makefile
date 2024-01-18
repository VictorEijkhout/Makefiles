################################################################
####
#### Makefile for Chaco installation
####
#### DOES NOT WORK: THE PETSC VERSION ONLY HAS A MAKEFILE
####
################################################################

PACKAGE = CHACO
URL = 
PACKAGEVERSION = git
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
.PHONY: default_install
default_install : configure install

GITREPO = https://bitbucket.org/petsc/pkg-chaco.git
include ${MAKEINCLUDES}/Make.git

