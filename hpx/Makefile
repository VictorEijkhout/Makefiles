################################################################
####
#### Makefile for HPX installation
####
################################################################

PACKAGE = HPX
URL = https://hpx.stellar-group.org/ \
    https://hpx-docs.stellar-group.org/branches/master/html/index.html
PACKAGEVERSION = git
MODE = omp

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
default_install : configure build public

GITREPO = https://github.com/STEllAR-GROUP/hpx.git
BRANCH = 1.10.0
include ${MAKEINCLUDES}/Make.git

