################################################################
####
#### Makefile for JULIA installation
####
#### build instructions:
#### https://github.com/JuliaLang/julia/blob/master/doc/src/devdocs/build/build.md
####
################################################################

PACKAGE = JULIA
URL = https://julialang.org/
PACKAGEVERSION = git
MODE = seq

include ${MAKEINCLUDES}/Make.info

## JULIA_INSTALL = ${STOCKYARD}/julia/installation-${TACC_SYSTEM}-${TACC_FAMILY_COMPILER}
info ::
	@echo "================ custom build"
	@echo "make default_install : all"
.PHONY: default_install
default_install :
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && ( \
	    echo "In-source build in $$srcdir" \
	     && rm -rf $$builddir && mkdir -p $$builddir \
	     && cd $$srcdir \
	     && echo " .. git reset and clean" \
	     && git reset --hard && git clean -x -f -d \
	     && echo prefix=$${prefixdir} > Make.user \
	     && echo " .. Make.user : " && echo "%%%%" && cat Make.user && echo "%%%%" \
	     && echo " .. starting make" \
	     && make --no-print-directory \
	     && make --no-print-directory install \
	    ) 2>&1 | tee $$installlog
##  BUILDROOT=$$builddir

.PHONY: projecttoml
build :: 
#projecttoml
projecttoml :
	@source ${MAKEINCLUDES}/names.sh \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && cd $$installdir \
	 && mkdir -p share/tacc \
	 && for line in \
	      "[extras]" \
	      "MPIPreferences = \"3da0fdf6-3ccc-4f1b-acd9-58baa6c99267\"" \
	      "" \
	      "[preferences.MPIPreferences]" \
	      "_format = \"1.0\"" \
	      "binary = \"system\"" \
	      "libmpi = \"libmpi\"" \
	      "version_string = \"Intel(R) MPI Library 2019 Update 9 for Linux* OS\n\"" \
	      "impl = \"IntelMPI\"" \
	      "version = v\"2019.9.0\"" \
	      "abi = \"MPICH\"" \
	      "mpiexec = \"ibrun\"" \
	    ; do \
	      echo $$line \
	    ; done > share/tacc/Project.toml

include ${MAKEINCLUDES}/Make.vars
HASBIN = 1
EXTRAINSTALLPATHS="JULIA_LOAD_PATH=share/tacc"
build ::
	@make --no-print-directory module

include ${MAKEINCLUDES}/Make.public
build ::
	@make public

GITREPO = https://github.com/JuliaLang/julia.git
BRANCH = v1.11.1
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

.PHONY: completely_clean
info ::
	@echo "make completely_clean : include .julia and such"
completely_clean : clean
	@/bin/rm -rf ${HOME}/.julia
	@source ${MAKEINCLUDES}/names.sh \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && rm -rf $$srcdir $$builddir $$installdir
