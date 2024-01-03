################################################################
####
#### Makefile for pdtoolkit
####
################################################################

PACKAGE = pdtoolkit
URL = https://www.pdtoolkit.org/
PACKAGEVERSION = 3.25.1
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

HASBIN = 1
NOLIB = 1
PREFIXOPTION = -prefix
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
default_install : configure build public

TGZURL = http://tau.uoregon.edu/pdt.tgz
include ${MAKEINCLUDES}/Make.download
