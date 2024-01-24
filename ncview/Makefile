################################################################
####
#### Makefile for NCVIEW installation
####
################################################################

PACKAGE = NCVIEW
URL = https://cirrus.ucsd.edu/ncview/
PACKAGEVERSION = 2.1.9
MODE = seq
MODULES = netcdf udunits
 
include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
HASBIN = 1
NOLIB = 1
CONFIGUREFLAGS = \
    --with-nc-config=${TACC_NETCDF_BIN}/nc-config \
    -with-udunits2_incdir=${TACC_UDUNITS_INC} \
    -with-udunits2_libdir=${TACC_UDUNITS_LIB}
default_install : 
	@export TACC_CC=$$( which ${TACC_CC} ) \
	 && make --no-print-directory configure build public

TGZURL = https://cirrus.ucsd.edu/~pierce/ncview/ncview-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
