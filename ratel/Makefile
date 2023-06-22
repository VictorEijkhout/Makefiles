################################################################
####
#### Makefile for Ratel installation
####
################################################################

URL = https://gitlab.com/micromorph/ratel
PACKAGE = RATEL
PACKAGEVERSION = git
CHMODR = 1
MODULES = petsc libceed
MODE = mpi

include ${MAKEINCLUDES}/Make.info

JCOUNT = 12
.PHONY: configure build
info ::
	echo "make configure build"
configure ::
	@echo "there is no separate configure rule"
build :: modules
	@( \
	source ${MAKEINCLUDES}/names.sh \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 && setcompilers && reportcompilers \
	 \
	 && echo "export CEED_DIR=[path to libCEED] PETSC_DIR=[path to PETSc] PETSC_ARCH=[PETSc arch] ENZYME_LIB=[path to Enzyme]" \
	 && export CEED_DIR=${TACC_LIBCEED_DIR} \
	 && export PETSC_DIR=${TACC_PETSC_DIR} && export PETSC_ARCH=${PETSC_ARCH} \
	 && echo "make with PETSC_DIR=$${PETSC_DIR}, PETSC_ARCH=$${PETSC_ARCH}, CEED_DIR=$${CEED_DIR}" \
	 && echo "into builddir: $${builddir}" \
	 && requirenonzero srcdir \
	 && cd $$srcdir \
	 && make --no-print-directory \
	        -B OPT=-g $${CCEXTRA} -j ${JCOUNT} V=1 \
	        BUILDDIR=$$builddir/ \
	        prefix=$$installdir \
	        lib install \
	 && ratelpc=$${installdir}/lib/pkgconfig/ratel.pc \
	 && if [ ! -f $$ratelpc ] ; then \
	        echo "was expecting to find ratel.pc at $$ratelpc" \
	    ; else \
	        echo ".. fixing ratel.pc" \
	         && [ 1 ] || sed -i \
	            -e '/prefix=/s?=.*?='$${builddir}'?' \
	            -e '/includedir=/s?=.*?='$${srcdir}/include'?' \
	            $${ratelpc} \
	         && cat $${ratelpc} \
	    ; fi \
	 && if [ "${CHMODR}" = "1" ] ; then \
	    echo "setting permissions open" \
	     && chmod -R g+rX,o+rX $${builddir} \
	; fi \
	 ) 2>&1 | tee install.log \
	 && echo "See install.log"

info ::
	@echo "make examples : copy and build examples"
.PHONY: examples
examples : petsc-installation
	@examples=examples-${PETSC_ARCH} \
	 && if [ "${TAU}" = "1" ] ; then \
	     examples=$${examples}-tau \
	    ; fi \
	 && rm -rf $${examples} \
	 && ( \
	    source ${MAKEINCLUDES}/names.sh \
	     && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} \
	     && ( cp -r $${srcdir}/examples $${examples} \
	         && cd $${examples} \
	         && builddir=$${builddir}-$${PETSC_ARCH} \
	         && if [ "${TAU}" = "1" ] ; then \
	                module load tau && export TAU_OPTIONS="$${TAU_OPTIONS} -optRevert" \
	                 && CCEXTRA="CC=tau_cc.sh" \
	                 && builddir=$${builddir}-tau \
	            ; fi \
	         && pcloc=$${installdir}/lib/pkgconfig \
	             && if [ -f $${pcloc}/ratel.pc ] ; then \
	                  echo "using $${pcloc}/ratel.pc:" \
	                   && cat $${pcloc}/ratel.pc \
	                ; else \
	                  echo "Can not find ratel.pc at: $${pcloc}" && exit 1 \
	                ; fi \
	         && export PKG_CONFIG_PATH=$${PKG_CONFIG_PATH}:$$pcloc \
	         && export CEED_DIR=$${PETSC_DIR}/$${PETSC_ARCH} \
	         && export RATEL_DIR=$${installdir} \
	         && echo "Make examples in `pwd`" \
	         && make --no-print-directory $${CCEXTRA} \
	        ) \
	    ) 2>&1 | tee $${examples}.log \
	 && echo && echo "Built examples in: $${examples}" \
	 && echo "see also $${examples}.log"

####
#### Run
####

SIZE = 44
.PHONY: submit
info :: 
	@echo "make submit [ SIZE=... from 44 3030 ] [ TAU= 0/1 ]"
submit :
	@echo "Submitting example run" \
	 && if [ "${TAU}" = "1" ] ; then \
	        echo ">>>>>" && echo Disabled: no tau examples && echo ">>>>" && exit 1 ; fi \
	 && ( \
	    source ${MAKEINCLUDES}/names.sh \
	     && MODE=mpi setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} \
	     && export examples_dir=`pwd`/examples-clx-intel-impi \
	     && if [ ! -d $${examples_dir} ] ; then \
	            echo "No such examples dir: $${examples_dir}" && exit 1 ; fi \
	     && echo "---- determining job file" \
	     && job=job${SIZE}.job \
	     && if [ ! -f $${job} ] ; then \
	            echo "Could not find job file: $${job}" && exit 1 \
	        ; else echo "---- Using job file: $${job}" ; fi \
	     && echo "---- submit" \
	     && LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$${installdir}/lib \
	        EXAMPLESDIR=$${examples_dir} sbatch $$job \
	    )

####
#### TAU
####

info ::
	@echo "make tau_process TAUDIR=..."
.PHONY: tau_process
TAUDIR = 
tau_process :
	@if [ ! -d "${TAUDIR}" ] ; then \
	    echo "TAUDIR not set or not found" && exit 1 ; fi \
	 && cd "${TAUDIR}" \
	 && rm -f tau.trc tau.edf align.trc align.edf \
	 && tau_treemerge.pl \
	 && tau_timecorrect tau.trc tau.edf align.trc align.edf \
	 && tau2slog2 align.trc align.edf -o yourprogram.slog2

####
#### Clean
####

include ${MAKEINCLUDES}/Make.clean
info ::
	@echo "================"
	@echo "make traceclean : remove logs and traces"
	@echo "make installclean : remove installations"
	@echo "make totalclean : everything clean"
.PHONY: installclean traceclean totalclean
installclean :: 
	@( \
	source ${MAKEINCLUDES}/names.sh \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} \
	 && rm -rf $${homedir}/build-* \
	 && rm -rf $${srcdir}/{build,lib} \
	)
errorclean :: clean
	@for j in $$(grep -l ERROR rateljob44.o50* | sed 's/rateljob44.o//' ) ; do \
	    echo "cleaning job $$j" \
	     && rm -rf tau_trace$${j} rateljob.o$${j} \
	 ; done
traceclean :: clean
	@rm -rf tau_trace*
clean ::
	@rm -f *~ *.log rateljob*.o*
totalclean :: traceclean installclean

####
#### Download stuff
####

include ${MAKEINCLUDES}/Make.git
GITREPO = https://gitlab.com/micromorph/ratel.git
