################################################################
####
#### Makefile for AUTOCONF installation
####
################################################################

PACKAGE = AUTOCONF
URL = https://ftp.gnu.org/gnu/autoconf/
PACKAGEVERSION = 2.72
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
NOLIB = 1
HASBIN = 1
include ${MAKEINCLUDES}/Make.install

.PHONY: default_install
default_install : configure build 

TGZURL = https://ftp.gnu.org/gnu/autoconf/autoconf-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean
