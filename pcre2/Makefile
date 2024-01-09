################################################################
####
#### Makefile for pcre2
####
#### https://www.pcre.org/
####
################################################################

PACKAGE = pcre2
URL = https://github.com/PCRE2Project/pcre2
PACKAGEVERSION = 10.42
MODE = seq

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install : configure build

TGZURL = https://github.com/PCRE2Project/pcre2/releases/download/pcre2-10.42/pcre2-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/PhilipHazel/pcre2.git
include ${MAKEINCLUDES}/Make.git
