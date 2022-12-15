################################################################
####
#### Makefile for Trilinos installation
####
################################################################

PACKAGE = trilinos
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 13-4-1
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODE = mpi
MODULES =
#boost swig

include ${MAKEINCLUDES}/Make.info

TRILINOS_VERSION=git
TRILINOS_HOME=${STOCKYARD}/trilinos
TRILINOS_SRC=${TRILINOS_HOME}/trilinos-${TRILINOS_VERSION}
TRILINOS_BUILD=${TRILINOS_HOME}/build
TRILINOS_INSTALL=${TRILINOS_HOME}/installation-${TACC_FAMILY_COMPILER}

info ::
	@echo "make configure install"
	@echo "    [ TRILINOS_VERSION=... (default=${TRILINOS_VERSION}) ]"
.PHONY: check_modules configure install config-help
check_modules :
	@if [ -z ${TACC_BOOST_DIR} ] ; then \
	    echo "Please load module: boost" && exit 1 ; fi
	@if [ -z ${TACC_SWIG_DIR} ] ; then \
	    echo "Please load module: swig" && exit 1 ; fi
configure :: check_modules
	@mkdir -p ${TRILINOS_HOME}
	@/bin/rm -rf ${TRILINOS_BUILD} && mkdir -p ${TRILINOS_BUILD}
	@export home=`pwd` \
	&& \
	( \
	cd ${TRILINOS_BUILD} \
	 && export PREFIXLOCATION==${TRILINOS_INSTALL} \
	 && export TRILINOS_LOCATION=${TRILINOS_HOME} && export VERSION=${TRILINOS_VERSION} \
	 && echo "Using install prefix=$${PREFIXLOCATION}" && echo \
	 && export COPTFLAGS="-O2" && export MKLFLAG="-mkl" \
	 && export PYTHON_VERSION=${TACC_FAMILY_PYTHON_VERSION} \
	 && export PYTHON_VERSION=$${PYTHON_VERSION%.*} \
	 && export PYTHON_LOAD_FLAG=${TACC_PYTHON_LIB}/libpython$${PYTHON_VERSION}m.so \
	 && export HAS_HDF5=0 && export HAS_NETCDF=0 && export HAS_MUMPS=0 \
	 && export HAS_SEACAS=0 && export HAS_STK=0 && export HAS_PYTHON=0 \
	 && export HAS_YAMLCPP=1 && source $${home}/../yaml-cpp/yamlcpp_vars.sh \
	 && source $${home}/trilinos_cmake.sh \
	 && ( \
	      echo "export TACC_TRILINOS_INC=$${INSTALL_PREFIX}/include" \
	       && echo "export TACC_TRILINOS_LIB=$${INSTALL_PREFIX}/lib64" \
	    ) > $$home/trilinos_vars.sh \
	 && echo "Set environment by sourcing: trilinos_vars.sh" \
	) 2>&1 | tee configure.log \

config-help :
	${TRILINOS_SRC}/configure --help

include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/trilinos/Trilinos/archive/refs/tags/trilinos-release-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = git@github.com:trilinos/Trilinos.git
include ${MAKEINCLUDES}/Make.git
