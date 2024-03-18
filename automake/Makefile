################################################################
####
#### Makefile for AUTOMAKE installation
#### https://automake.org/index.html
####
################################################################

PACKAGE = AUTOMAKE
URL = https://ftp.gnu.org/gnu/automake/
PACKAGEVERSION = 1.16.5
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
NOLIB = 1
HASBIN = 1
include ${MAKEINCLUDES}/Make.install

.PHONY: default_install
default_install : configure build 

TGZURL = https://ftp.gnu.org/gnu/automake/automake-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean
