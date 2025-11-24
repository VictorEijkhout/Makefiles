################################################################
####
#### Makefile for LIBTOOL installation
####
################################################################

PACKAGE = LIBTOOL
URL = https://ftp.gnu.org/gnu/libtool/
PACKAGEVERSION = 2.4.7
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
NOLIB = 1
HASBIN = 1
include ${MAKEINCLUDES}/Make.install

.PHONY: default_install
default_install : configure build 

DOWNLOADURL = https://ftp.gnu.org/gnu/libtool/libtool-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean
