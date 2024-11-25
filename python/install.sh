#!/bin/bash

## https://devguide.python.org/getting-started/setup-building/#setup
## https://github.com/python/cpython/issues/122412

function usage () {
    echo "Usage: $0 [ -h ] [ -v 3.12345 (default: ${pythonver}) ]"
    echo "    [ -5 (hdf5 only) ] [ -m (mpi only) ] [ -n (numpy only) ] "
    echo "    [ -o (others only: ${others}) ]"
    echo "    [ --prefix prefixdir ]"
    echo "    no options: install everything"
}

installhdf=1
installmpi=1
installnumpy=1
installpython=1
installothers=1
others="paramiko setuptools"
prefixdir=

pythonver=3.12.4
#3.11.0
#

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
    elif [ $1 = "-v" ] ; then
	shift && pythonver=$1 && shift
    else
	echo "Unknown option: <<$1>>" && exit 1
    fi
done

if [ -z "${TACC_HDF5_DIR}" ] ; then
    echo "Needs hdf5 module loaded" && exit 1
fi

pymacrover=${pythonver%%.*}
pyminiver=${pythonver#*.} && pyminiver=${pyminiver%.*}
pymicrover=${pythonver##*.}
echo && echo "Installing ${pymacrover}.${pyminiver}.${pymicrover}" && echo

pythondir=${STOCKYARD}/python
if [ -z "${prefixdir}" ] ; then 
    prefixdir=${pythondir}/installation-${pythonver}-${TACC_SYSTEM}-${TACC_FAMILY_COMPILER}-${TACC_FAMILY_COMPILER_VERSION}
fi
pkgprefix=${prefixdir}/lib/python${pymacrover}.${pyminiver}/site-packages/

if [ ! -z "${installpython}" ] ; then 
    rm -rf ${HOME}/{.local,.cache}
    echo && echo "Building python ${pythonver}" && echo
    cd ${pythondir}
    # if [ ! -f Python-${pythonver}.tgz ] ; then 
    # 	echo && echo "first downloading tgz" && echo
    # 	wget https://www.python.org/ftp/python/${pythonver}/Python-${pythonver}.tgz
    # fi
    # if [ ! -f Python-${pythonver} ] ; then 
    # 	echo && echo "Untar" && echo
    # 	tar fxz Python-${pythonver}.tgz
    # fi
    # echo && echo "Remove previous installation" && echo
    # rm -rf ${prefixdir}

    cd ${pythondir}/Python-${pythonver}

    module load sqlite || exit 1
    export CC=${TACC_CC} && export CXX=${TACC_CXX}
    if [ "${TACC_COMPILER_FAMILY}"  = "intel" ] ; then 
	export LDFLAGS="-DLDFLAGS -L${TACC_MKL_LIB:?MISSING_MKL_LIB} -L${TACC_INTEL_LIB:?MISSING_INTEL_LIB}"
    else 
	echo "find /opt/intel/oneapi/compiler/2024.0/lib/libintlc.so.5"
	export LDFLAGS="-DLDFLAGS \
-L${TACC_MKL_LIB:?MISSING_MKL_LIB}       -Wl,-rpath=${TACC_MKL_LIB} \
-L/opt/intel/oneapi/compiler/2024.0/lib  -Wl,-rpath=/opt/intel/oneapi/compiler/2024.0/lib \
-lintlc"
    fi
    echo && echo "Configuring" && echo
    ##./configure --help 
    ./configure --prefix=${prefixdir} \
		--enable-optimizations \
		--with-ensurepip=install \
		2>&1 | tee ${pythondir}/configure.log
    echo && echo "Making" && echo
    ( make -j 1 && echo && echo "Make install" && echo && make -j 1 install ) 2>&1
    ## | tee ${pythondir}/install.log
fi

##
## set paths for the rest of the installation
##
if [ ! -d "${pkgprefix}" ] ; then 
    echo "ERROR: not finding site-packages: ${pkgprefix}" && exit 1 
fi
cat >${pythondir}/load.sh <<EOF
export TACC_PYTHON_DIR=${prefixdir}
export PATH=${prefixdir}/bin:/home1/00434/eijkhout/.local/bin:${PATH}
export PYTHONPATH=${pkgprefix}:${PYTHONPATH}
EOF
source ${pythondir}/load.sh
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
## https://docs.scipy.org/doc/scipy/building/blas_lapack.html
##
numpyver=2.0.1
## 1.21.0
numpygit=${pythondir}/numpy-git
numpydir=${pythondir}/numpy-git
scipyver=1.14.0
scipygit=${pythondir}/scipy-git

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
    # for p in $( echo ${PKG_CONFIG_PATH} | tr ':' ' ' ) ; do
    # 	echo "PKG_CONFIG_PATH contains: $p"
    # 	ls ${p}/*.pc
    # done
    # for p in $( echo ${CMAKE_PREFIX_PATH} | tr ':' ' ' ) ; do
    # 	echo "CMAKE_PREFIX_PATH contains $p"
    # 	find $p -name \*.cmake
    # done
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

exit 0

################################################################
##
## problem with Intel:
##

Traceback (most recent call last):
  File "/work2/00434/eijkhout/stampede3/python/Python-3.12.4/Lib/multiprocessing/process.py", line 314, in _bootstrap
    self.run()
  File "/work2/00434/eijkhout/stampede3/python/Python-3.12.4/Lib/multiprocessing/process.py", line 108, in run
    self._target(*self._args, **self._kwargs)
  File "/work2/00434/eijkhout/stampede3/python/installation-3.12.4-intel/lib/python3.12/concurrent/futures/process.py", line 251, in _process_worker
    call_item = call_queue.get(block=True)
                ^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/work2/00434/eijkhout/stampede3/python/Python-3.12.4/Lib/multiprocessing/queues.py", line 102, in get
    with self._rlock:
  File "/work2/00434/eijkhout/stampede3/python/Python-3.12.4/Lib/multiprocessing/synchronize.py", line 95, in __enter__
    return self._semlock.__enter__()
           ^^^^^^^^^^^^^^^^^^^^^^^^^
