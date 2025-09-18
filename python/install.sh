#!/bin/bash

## https://devguide.python.org/getting-started/setup-building/#setup
## https://github.com/python/cpython/issues/122412

function usage () {
    echo "Usage: $0 [ -h ] [ -v 3.12345 (default: ${pythonver}) ]"
    echo "    [ -5 (hdf5 only) ] [ -m (mpi only) ] [ -n (numpy only) ] "
    echo "    [ -o (others only: ${others}) ]"
    echo "    [ --build builddir ]"
    echo "    [ --prefix prefixdir ]"
    echo "    [ --srcdir srcdir ]"
    echo "    no options: install everything"
}

installhdf=1
installmpi=1
installnumpy=1
installpython=1
installothers=1
others="paramiko setuptools"
builddir=
prefixdir=
srcdir=

pythonver=3.12.5

while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 0
    elif [ $1 = "-5" ] ; then
	shift && installhdf=1 && installpython= && installmpi= && installnumpy= && installothers=
    elif [ $1 = "-m" ] ; then
	shift && installhdf= && installpython= && installmpi=1 && installnumpy= && installothers=
    elif [ $1 = "-n" ] ; then
	shift && installhdf= && installpython= && installmpi= && installnumpy=1 && installothers=
    elif [ $1 = "-o" ] ; then
	shift && installhdf= && installpython= && installmpi= && installnumpy= && installothers=1
    elif [ $1 = "--prefix" ] ; then
	shift && prefixdir=$1 && shift
    elif [ $1 = "--build" ] ; then
	shift && builddir=$1 && shift
    elif [ $1 = "--srcdir" ] ; then
	shift && srcdir=$1 && shift
    elif [ $1 = "-v" ] ; then
	shift && pythonver=$1 && shift
    else
	echo "Unknown option: <<$1>>" && exit 1
    fi
done
if [ -z "${srcdir}" ] ; then 
    srcdir=${STOCKYARD}/python/python-${pythonver}
    if [ ! -d "${srcdir}" ]  ; then 
	echo "Need to download python-${pythonver} first!" && exit 1
    fi
fi
if [ -z "${TACC_HDF5_DIR}" ] ; then
    echo "Needs hdf5 module loaded" && exit 1
fi

##
## Version handing
##
pymacrover=${pythonver%%.*}
pyminiver=${pythonver#*.} && pyminiver=${pyminiver%.*}
pymicrover=${pythonver##*.}
echo && echo "Installing ${pymacrover}.${pyminiver}.${pymicrover}" && echo

##
## Prefix
##
if [ -z "${prefixdir}" ] ; then 
    # standard pythondir for local install only
    pythondir=${STOCKYARD}/python
    prefixdir=${pythondir}/installation-${pythonver}-${TACC_SYSTEM}-${TACC_FAMILY_COMPILER}-${TACC_FAMILY_COMPILER_VERSION}
    mkdir -p "${prefixdir}"
fi
pkgprefix=${prefixdir}/lib/python${pymacrover}.${pyminiver}/site-packages/
mkdir -p "${pkgprefix}"
echo "${pkgprefix}" > ${pkgprefix}/tacc.pth

##
## Configure and install python
##
if [ ! -z "${installpython}" ] ; then 

    echo && echo "================ Install basic python ================" && echo

    if [ -z "${srcdir}" ] ; then echo "Zero variable srcdir" ; exit 1 ; fi
    echo "Using srcdir=${srcdir}"
    cd ${srcdir}

    export CC=${TACC_CC} && export CXX=${TACC_CXX}
    if [ "${TACC_COMPILER_FAMILY}"  = "intel" ] ; then 
	export LDFLAGS="-DLDFLAGS -L${TACC_MKL_LIB:?MISSING_MKL_LIB} -L${TACC_INTEL_LIB:?MISSING_INTEL_LIB}"
    else 
	echo "find libintlc.so"
	# /opt/intel/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64_lin/libintlc.so
	if [ -z "${TACC_MKL_DIR}" ] ; then echo "Please load MKL" && exit 1 ; fi
	libintlc=$( find ${TACC_MKL_DIR}/../.. -name libintlc.so 2>/dev/null | cut -d ' ' -f 1 )
	echo "found lib: ${libintlc}"
	libintlc=${libintlc%%libintlc.so*}
	echo "found dir: ${libintlc}"
	export LDFLAGS="-DLDFLAGS \
-L${TACC_MKL_LIB:?MISSING_MKL_LIB} -Wl,-rpath=${TACC_MKL_LIB} \
-L${libintlc}                      -Wl,-rpath=${libintlc}     \
-lintlc"
    fi
    echo && echo "Configuring" && echo

    ./configure --prefix=${prefixdir} \
		--disable-test-modules \
		--enable-optimizations \
		--with-ensurepip=install \
	--with-openssl=/usr --with-openssl-rpath=/usr/lib64/openssl

    echo && echo "Making" && echo

    echo " .. WARNING weird libgdbm patch"
    rm -rf $(pwd)/lib64 && mkdir -p lib64 && ln -s /usr/lib64/libgdbm.so.6 ./lib64/libgdbm.so
    export LD_LIBRARY_PATH=$(pwd)/lib64:${LD_LIBRARY_PATH}

    echo " .. parallel make"
    make -j 12 LDFLAGS_NODIST=-L$( pwd )/lib64
    echo && echo "Make install" && echo
    make -j 5 install
fi

##
## set paths for the rest of the installation
##
if [ ! -d "${pkgprefix}" ] ; then 
    echo "ERROR: not finding site-packages: ${pkgprefix}" && exit 1 
fi
cat >${builddir:?MISSING_BUILD_DIR}/load.sh <<EOF
export TACC_PYTHON_DIR=${prefixdir}
export PATH=${prefixdir}/bin:/home1/00434/eijkhout/.local/bin:${PATH}
export PYTHONPATH=${pkgprefix}:${PYTHONPATH}
EOF
source ${builddir}/load.sh
echo "Now using python: $( which python3 ), deduced version: $( python3 --version )"
export XDG_CACHE_HOME=$( pwd )/pipcache

##
## new pip
##
echo && echo "Upgrading pip" && echo
python3 -m pip install --upgrade pip

##
## numpy & scipy
## https://github.com/numpy/numpy
## https://github.com/scipy/scipy.git
## https://docs.scipy.org/doc/scipy/building/blas_lapack.html
##
numpyver=2.0.1
## 1.21.0
numpygit=${builddir}/numpy-git
if [ ! -d ${numpygit} ] ; then
    ( cd ${builddir} && git clone https://github.com/numpy/numpy numpy-git )
fi
scipyver=1.14.0
scipygit=${builddir}/scipy-git
if [ ! -d ${scipygit} ] ; then
    ( cd ${builddir} && git clone https://github.com/scipy/scipy.git scipy-git )
fi

if [ ! -z "${installnumpy}" ] ; then
    rm -rf ${pkgprefix}/numpy*
    rm -rf ${HOME}/.cache/pip ${XDG_CACHE_HOME}

    ##
    ## first numpy install
    ##
    if [ ! -z "" ] ; then
	echo && echo "Building Numpy ${numpyver}" && echo
	pushd ${numpygit}
	git fetch --all --tags --prune
	git checkout tags/v${numpyver}
	## https://github.com/numpy/numpy/pull/26802
	pip3 install . --target=${pkgprefix} \
             -Csetup-args=-Dblas=${blas} -Csetup-args=-Dlapack=${blas}
	## python3 -m pip install --target=${pkgprefix} .
	## \
	    ## -Csetup-args=-Dblas-order=mkl,blis,openblas \
	    ## -Csetup-args=-Dallow-noblas=false \
	    ## -Csetup-args=-Dlapack-order=mkl,openblas,lapack
	popd
    fi
    
    ##
    ## now scipy
    ##
    echo && echo "Building Scipy ${scipyver}" && echo
    pushd ${scipygit}
    git submodule update --init
    git fetch --all --tags --prune
    git checkout tags/v${scipyver}
    pip3 --version
    blas=mkl-sdl
    #blas=blis
    pip3 install . --target=${pkgprefix} \
	 -Csetup-args=-Dblas=${blas} -Csetup-args=-Dlapack=${blas}
    ## -Csetup-args=-Dblas=mkl -Csetup-args=-Dlapack=mkl
	 # -Csetup-args=-Dblas-order=mkl,blis,openblas \
	 # -Csetup-args=-Dallow-noblas=false \
	 # -Csetup-args=-Dlapack-order=mkl,openblas,lapack
    popd
    # pip3 install --target=${pkgprefix} \
    # 	    scipy \
    #     -Csetup-args=-Dblas=mkl -Csetup-args=-Dlapack=mkl
fi

##
## hdf5
##
if [ ! -z "${installhdf}" ] ; then
    echo && echo "Building h5py" && echo
    HDF5_DIR=${TACC_HDF5_DIR} \
    HDF5_VERSION=${TACC_HDF5_VERSION} \
	pip3 install --target=${pkgprefix} --no-binary=h5py h5py
    ## CC="mpicc" HDF5_MPI="ON" HDF5_DIR=/path/to/parallel-hdf5 pip3 install --target=${pkgprefix} --no-binary=h5py h5py
fi

##
## mpi4py
##
if [ ! -z "${installmpi}" ] ; then
    echo && echo "Building mpi4py" && echo
    python3 -m pip install --target=${pkgprefix} \
	    mpi4py
fi

##
## others
##
function install () {
    package=$1
    echo && echo "Building ${package}" && echo
    python3 -m pip install --target=${pkgprefix} \
	    ${package}
}

if [ ! -z "${installothers}" ] ; then
    for o in ${others} ; do
	install $o
    done
fi

chmod -R g+rX,o+rX ${prefixdir}

moduledir=${WORK}/modulefiles/Compiler/${TACC_FAMILY_COMPILER}/${TACC_FAMILY_COMPILER_VERSION}/python3
mkdir -p ${moduledir}
cat >${moduledir}/${pythonver}.lua <<EOF
setenv( "TACC_PYTHON_DIR",  "${prefixdir}" )
setenv( "TACC_PYTHON3_DIR", "${prefixdir}" )

prepend_path( "PATH",       "${prefixdir}/bin" )
prepend_path( "PYTHONPATH", "${pkgprefix}" )
EOF

echo && echo "================ Installation finished"
echo "in: ${prefixdir}"
ls -ld ${prefixdir}

