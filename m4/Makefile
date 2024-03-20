################################################################
####
#### Makefile for M4 installation
####
################################################################

PACKAGE = M4
URL = https://ftp.gnu.org/gnu/m4/
PACKAGEVERSION = 1.4.19
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
NOLIB = 1
HASBIN = 1
include ${MAKEINCLUDES}/Make.install

.PHONY: default_install
default_install : configure build 

TGZURL = https://ftp.gnu.org/gnu/m4/m4-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean
