info ::
	echo "Available targets:"

.PHONY: test
info ::
	@echo "make test [STEP=12345 (default=${STEP})]"
STEP = 36
test :
	export srcdir=`pwd`/step-${STEP} \
	&& export builddir=`pwd`/build \
	&& export installdir=`pwd`/install \
	&& rm -rf $${builddir} $${installdir} \
	&& mkdir -p $${builddir} \
	&& export prmfile=step-${STEP}.prm \
	&& if [ -f $${prmfile} ] ; then \
	     cp $${prmfile} $${builddir} \
	     ; fi \
	&& mkdir -p $${installdir} \
	&& ( \
	     cd $${builddir} \
	     && cmake \
	        -D CMAKE_INSTALL_PREFIX=$${installdir} \
	        -D DEAL_II_DIR=${TACC_DEALII_DIR} \
	        $${srcdir} \
	     && make \
	   )
