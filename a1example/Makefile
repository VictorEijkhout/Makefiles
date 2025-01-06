################################################################
####
#### Example makefile
####
################################################################

##
## PACKAGE name can be in any case:
## we normalize to all lowercase anyway
PACKAGE = siesta

##
## URL goes into the modulefile,
## display with "make url"
## or open with "make homepage"
##
URL = https://siesta-project.org/siesta/index.html

##
## PACKAGEVERSION as it goes into the modulefile
## sometimes you need to declare variants such as
## PACKAGEVERSIONDOWNLOAD for the tgz line below
##
PACKAGEVERSION = 5.2.1

##
## MODE determines where in the module hierarchy the package is placed:
## - core : no compiler dependece
## - omp / seq : only dependent on compiler, not on MPI
## - mpi : dependent on MPI
##
MODE = mpi

##
## MODULES is a list of prerequisite modules for building
## - mkl : only loaded for gcc; for nvidia this loads nvpl, for intel it's a no-op
## - these modules do not go into the "depends_on" clause of the module file
##
MODULES = fftw3 mkl 
##

##
## The first of the include files.
## this one is absolutely required
##
include ${MAKEINCLUDES}/Make.info

####
####
#### Configure and build
####
####

##
## CMAKE installation needs the following include files:
##
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

##
## CMAKEFLAGS and related are for the cmake-based build
## Seeral obvious cmake flags are provided by default
##
CMAKEFLAGS = \
    -D PACKAGE_USE_OTHERPACKAGE=ON
##
## cmake build type is RelWithDebug by default
##
CMAKEBUILDDEBUG = 1

##
## if CMakeLists.txt is not at root level:
##
CMAKESUBDIR=src

##
## CMake has its own logic whether to generate a `lib' or `lib64' directory.
## Since this goes into an environment variable you don't care.
## However, if you hard rely on `lib' being there, set:
##
LINKLIB64toLIB = 1

##
## Autotools uses the lines
## !! incompatible with cmake: use one or the other !!
##
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install

##
## if configure / make files are not at root level:
##
CONFIGURESUBDIR=src
MAKESUBDIR=src

##
## The scripts are pretty clever about `autoreconf' and such,
## however you can always use:
##
BEFORECONFIGURECMDS = autoreconf -i

##
## Customization. any non-zero value counts as true; default is false
## - HASBIN : is there a bin directory? 
## - NOLIB : is there no lib directory?
## - NOINCLUDE : is there no include directory? 
HASBIN = 1
NOLIB =
NOINCLUDE =

##
## Customization settings:
## - CMAKEPREFIXPATHSET : nonzero causes the CMAKE_PREFIX_PATH to be set
## - PKGCONFIG : path for ".pc" files relative to installation root
## - PKGCONFIGLIB : path for ".pc" files relative to lib/lib64 dir
## - CPTOINSTALLDIR : extra files/directories to copy to the installation
## - INCLUDELOC : non-standard name/location of the include directory
##
PKGCONFIGLIB = pkgconfig
CPTOINSTALLDIR = Examples

##
## Convention: have a "default_install" rule that does obvious stuff.
## For some packages (like hdf5) you want to have separate "seq/par" rules
##
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

##
## parallel make: set
## JCOUNT=24
##

####
####
#### Download / clone
####
####

##
## Download of versioned releases
## Specify a TGZURL path
## - this can end in tgz or tar.gz; it gets renamed to a standard name anyway
## - use PACKAGEVERSION here as much as possible
## - for zip files use ZIPURL
##
TGZURL = https://gitlab.com/pkg-project/pkg/-/archive/${PACKAGEVERSION}/pkg-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

##
## Clone of git repo
## Specify a GITREPO path
## This rule can coexist with the TGZURL download above,
## however, you then need to set the package version explicitly:
##   make clone PACKAGEVERSION=git
## - SUBMODULE : set to nonzero to have submodules loaded
##
GITREPO = https://gitlab.com/pkg-project/pkg.git
SUBMODULE = 1
include ${MAKEINCLUDES}/Make.git
