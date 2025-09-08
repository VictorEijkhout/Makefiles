################################################################
####
#### Makefile for Ncl installation
####
################################################################

PACKAGE = ncl
ABOUT = NCAR Command Language
PACKAGEVERSION = 6.6.2
PACKAGEVERSION_NODOT = 662
URL = https://www.ncl.ucar.edu/
DOCURL = https://www.ncl.ucar.edu/Download/build_from_src.shtml
MODE = seq

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "make system_test : can you build ncl?"
.PHONY: system_test
system_test :
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && requirenonzero prefixdir && mkdir -p $${prefixdir} && rm -rf $${prefixdir}/* \
	 && requirenonzero builddir && rm -rf $$builddir && mkdir -p $$builddir \
	 && requirenonzero logfilesdir \
	 \
	 && cd $${srcdir}/config &&  make -f Makefile.ini \
	 && ./ymake -config `pwd` 2>&1 | tee parse \
	 && if [ $$( cat parse | wc -l ) -gt 0 ] ; then \
	        echo "System not recognized" && cat parse \
	    ; else \
	        echo "System recognized:" \
	         && grep SYSTEM_INCLUDE Makefile \
	    ; fi
.PHONY: configure build
configure : system_test
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && requirenonzero prefixdir && mkdir -p $${prefixdir} && rm -rf $${prefixdir}/* \
	 && requirenonzero builddir && rm -rf $$builddir && mkdir -p $$builddir \
	 && requirenonzero logfilesdir \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 \
	 && ( \
	    cat Makefile.template \
	            | sed -e 's?theinstalldir?'$${prefixdir}'?' \
	            > $${srcdir}/common/Makefile \
	     && cd $${srcdir}/common \
	     && echo && echo "==== Using configuration:" && make Info \
	    ) 2>&1 | tee $$configurelog
build :
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && requirenonzero prefixdir && mkdir -p $${prefixdir} && rm -rf $${prefixdir}/* \
	 && requirenonzero builddir && rm -rf $$builddir && mkdir -p $$builddir \
	 && requirenonzero logfilesdir \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 \
	 && mkdir -p $${prefixdir} \
	 && ( \
	    cd $${srcdir}/common \
	     && echo && echo "==== Starting installation into $${prefixdir} with:" \
	     && make Info \
	          INSTALLDIR=$${prefixdir} \
	     && echo \
	     && make Makefiles clean includes depend all install \
	          INSTALLDIR=$${prefixdir} \
	    ) 2>&1 | tee $${installlog} \
	 && echo && echo "Installation finished in: $$prefixdir" \
	 && echo "(see: $$installlog)" && echo

NOLIB = 1
HASBIN = 1
include ${MAKEINCLUDES}/Make.vars

info :: 
	@echo "make default_install : configure build module"
.PHONY: default_install
default_install : configure build module

## include ${MAKEINCLUDES}/Make.configure
## include ${MAKEINCLUDES}/Make.install

TGZURL = https://github.com/NCAR/ncl/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/NCAR/ncl.git
include ${MAKEINCLUDES}/Make.git

