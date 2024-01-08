################################################################
####
#### Makefile for HPX installation
####
################################################################

PACKAGE = HPX
URL = https://hpx.stellar-group.org/ \
    https://github.com/STEllAR-GROUP/hpx \
    https://hpx-docs.stellar-group.org/branches/master/html/index.html
PACKAGEVERSION = 1.9.1
MODE = seq
MODULES = boost

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
CMAKEFLAGS = \
    -D Boost_ROOT=${TACC_BOOST_DIR} \
    -D HPX_WITH_MALLOC=system \
    -D HPX_WITH_FETCH_ASIO=ON
HASBIN = 1
PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.cbuild
default_install : configure build public

TGZURL = https://github.com/STEllAR-GROUP/hpx/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/STEllAR-GROUP/hpx.git
# BRANCH = 1.10.0
include ${MAKEINCLUDES}/Make.git

