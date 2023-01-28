################################################################
####
#### Makefile for Conquest installation
####
#### https://conquest.geodynamics.org/
####
################################################################

PACKAGE = CONQUEST
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = GIT
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODE = mpi
MODULES = fftw3 

include ${MAKEINCLUDES}/Make.info

build :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} "" ${PACKAGEBASENAME} \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 && setmpicompilers \
	    \
	 && rm -rf $$builddir && mkdir -p $$builddir \
	 && ( \
	    cd $${builddir} \
	     && reportnames && echo \
	     && reportcompilers && echo \
	    ) 2>&1 | tee $$configurelog \
	 && make --no-print-directory varsmodule


GITREPO = git clone https://github.com/OrderN/CONQUEST-release
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

