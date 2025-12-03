################################################################
####
#### Makefile for HMMER installation
####
################################################################

PACKAGE = HMMER
ABOUT = 
URL = 
PACKAGEVERSION = 3.4
MODE = core
MODULES = 

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
NOLIB = 1
NOINC = 1
HASBIN = 1

.PHONY: default_install
default_install : configure build 

DOWNLOADURL = http://eddylab.org/software/hmmer/hmmer.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean
