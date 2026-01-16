################################################################
####
#### Makefile for pcre2
####
#### https://www.pcre.org/
####
################################################################

PACKAGE = pcre2
ABOUT = Perl-compatible Regular Expressions
URL = https://github.com/PCRE2Project/pcre2
PACKAGEVERSION = 10.46
MODE = seq

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET=1
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install : configure build

DOWNLOADURL = https://github.com/PCRE2Project/pcre2/archive/refs/tags/pcre2-${PACKAGEVERSION}.tar.gz
## https://github.com/PCRE2Project/pcre2/releases/download/pcre2-10.42/pcre2-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/PhilipHazel/pcre2.git
include ${MAKEINCLUDES}/Make.git
