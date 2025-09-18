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
	        echo "System not recognized" && cat parse && exit 1 \
	    ; else \
	        export system=$$( grep SYSTEM_INCLUDE Makefile | cut -d '=' -f 2 | tr -d '\"' | tr -d " " ) \
	         && echo "System recognized: <<$${system}>>" \
	    ; fi \
	 && requirenonzero scriptdir \
	 && cp $${scriptdir}/config_$${system}_${TACC_FAMILY_COMPILER} $${srcdir}/config/$${system}
.PHONY: makefile configure build
makefile :
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && requirenonzero prefixdir \
	 && requirenonzero builddir \
	 && requirenonzero scriptdir \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 \
	 && cd $${srcdir} \
	 && ( echo "Creating makefile" \
	     && ( echo "hm?" >/dev/null \
	        && echo "overwrite"  >/dev/null && echo "y" \
	        && echo "ack"        >/dev/null && echo \
	        && echo "ack"        >/dev/null && echo \
	        && echo "build ncl"  >/dev/null && echo "y" \
	        && echo "prefix"     >/dev/null && echo "$${prefixdir}" \
	        && echo "tmp space"  >/dev/null && echo "$${builddir}" \
	        && echo "ncdf4"      >/dev/null && echo "n" \
	        && echo "hdf4"       >/dev/null && echo "n" && echo "n" \
	        && echo "triangle"   >/dev/null && echo "n" \
	        && echo " .. ncf4"   >/dev/null && echo "n" && echo "n" \
	        && echo "gdal"       >/dev/null && echo "y" \
	        && echo "eemd"       >/dev/null && echo "n" \
	        && echo "udu2"       >/dev/null && echo "y" \
	        && echo "vis5"       >/dev/null && echo "n" \
	        && echo "eos2"       >/dev/null && echo "n" \
	        && echo "hdf5"       >/dev/null && echo "y" && echo "n" \
	        && echo "grib"       >/dev/null && echo "n" \
	        && echo "lib search" >/dev/null && echo "/usr/lib64/X11/" \
	        && echo " .. inc "   >/dev/null && echo "/usr/include/X11" \
	        && echo "back/save"  >/dev/null && echo "n" && echo "y" \
	        ) | ./Configure -v \
	     && make Info \
	    ) 2>&1 | tee $${configurelog} \
	 && echo && echo "see $${configurelog}" && echo
everything :
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && requirenonzero prefixdir && mkdir -p $${prefixdir} && rm -rf $${prefixdir}/* \
	 && requirenonzero builddir && rm -rf $$builddir && mkdir -p $$builddir \
	 && requirenonzero installlog \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 \
	 && mkdir -p $${prefixdir} \
	 && ( \
	    cd $${srcdir} && make Everything \
	    ) 2>&1 | tee $${installlog}

NOLIB = 1
HASBIN = 1
include ${MAKEINCLUDES}/Make.vars

info :: 
	@echo "make default_install : configure build module"
.PHONY: default_install
default_install : system_test makefile everything varsmodule public

TGZURL = https://github.com/NCAR/ncl/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/NCAR/ncl.git
include ${MAKEINCLUDES}/Make.git

################################################################
################################################################
################################################################
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
