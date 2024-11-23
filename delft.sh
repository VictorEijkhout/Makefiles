#!/bin/bash

function load_or_install_module () {
    m=$1
    echo "testing $m"
    module is-avail $m
    if [ $? -eq 1 ] ; then 
	echo " .. installing $m"
	pushd $m 
	case $m in 
	    ( hdf5 ) make seq ;;
	    ( netcdf ) make seq ;;
	    ( petsc ) ./install_small.sh ;;
	    ( * ) make default_install ;;
	esac
	popd
    fi
    module load $m
    # module -t list
}

prereqs="pcre2 bison swig zlib"
delft="hdf5 netcdf jsonc sqlite proj gdal metis"
petsc="eigen fftw3 petsc"
for m in ${prereqs} ${delft} ${petsc} ; do
    load_or_install_module $m
done

cmdline="module -t avail ${delft} petsc"
echo "What is available?"
echo $cmdline
eval $cmdline
