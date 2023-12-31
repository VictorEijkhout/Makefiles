################################################################
####
#### Makefile for HPX installation
####
################################################################

PACKAGE = HPX
URL = https://hpx.stellar-group.org/ \
    https://hpx-docs.stellar-group.org/branches/master/html/index.html
PACKAGEVERSION = git
MODE = seq
MODULES = boost

include ${MAKEINCLUDES}/Make.info

# PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.cmake
CMAKEFLAGS = \
    -D Boost_ROOT=${TACC_BOOST_DIR} \
    -D HPX_WITH_MALLOC=system \
    -D HPX_WITH_FETCH_ASIO=ON
include ${MAKEINCLUDES}/Make.cbuild
default_install : configure build public

GITREPO = https://github.com/STEllAR-GROUP/hpx.git
# BRANCH = 1.10.0
include ${MAKEINCLUDES}/Make.git

