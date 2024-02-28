################################################################
####
#### Makefile for Parallel-Netcdf installation
####
################################################################

PACKAGE = PARALLELNETCDF
URL = https://parallel-netcdf.github.io/ \
    https://github.com/Parallel-NetCDF/PnetCDF
PACKAGEVERSION = git
# 1.12.3
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

BEFORECONFIGURECMDS = autoreconf -i
include ${MAKEINCLUDES}/Make.configure

HASBIN=1
# CMAKEPREFIXPATHSET = 1
PKGCONFIGLIB = pkgconfig

##
## NOTE the enable_pnetcdf flag is for parallel-netcdf 1.9
##
CMAKEFLAGS = 
include ${MAKEINCLUDES}/Make.install
.PHONY: default_install
info ::
	@echo "make default_install: parallel version, also available: seq"
default_install : configure build

TGZURL = https://github.com/Parallel-NetCDF/PnetCDF/archive/refs/tags/checkpoint.${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/Parallel-NetCDF/PnetCDF.git
include ${MAKEINCLUDES}/Make.git
include ${MAKEINCLUDES}/Make.clean
