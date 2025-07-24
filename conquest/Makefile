################################################################
####
#### Makefile for Conquest installation
####
#### https://conquest.geodynamics.org/
####
################################################################

PACKAGE = CONQUEST
PACKAGEVERSION = git
MODE = mpi
MODULES = fftw3 
# xclib

include ${MAKEINCLUDES}/Make.info

.PHONY: configure build
info ::
	@echo "make build (there is no configure stage)"
	@echo "    [ NOCLEAN= nonzero for no make clean ]"
configure :
	@echo && echo "there is no configure rule" && echo
build :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} "" ${PACKAGEBASENAME} \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 && setmpicompilers \
	    \
	 && rm -rf $$installdir \
	 && cp system.make.${TACC_SYSTEM} $$srcdir/src/system.make \
	 && cp system.make.${TACC_SYSTEM} $$srcdir/tools/BasisGeneration/system.make \
	 && ( \
	    cd $${srcdir}/src \
	     && reportnames && echo \
	     && reportcompilers && echo \
	     && if [ ! -z "${NOCLEAN}" ] ; then \
	            make very_clean && ( cd ../tools/BasisGeneration && make clean ) \
	        ; fi \
	     && make \
	     && mkdir -p $$installdir/bin \
	     && cp Conquest $$installdir/bin \
	     \
	     && ( cd ../tools/BasisGeneration && make ) \
	     && cp ../bin/* $$installdir/bin \
	    ) 2>&1 | tee $$installlog \
	 && make --no-print-directory varsmodule
include ${MAKEINCLUDES}/Make.vars
include ${MAKEINCLUDES}/Make.public

GITREPO = https://github.com/OrderN/CONQUEST-release
include ${MAKEINCLUDES}/Make.git

ZIPURL = https://github.com/OrderN/CONQUEST-release/archive/master.zip
TGZRENAME = mv master.zip conquest-git.zip
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean

