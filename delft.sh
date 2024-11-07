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
	    ( * ) make default_install ;;
	esac
	popd
    fi
    module load $m
    # module -t list
}

for m in \
        pcre2 bison swig \
        zlib hdf5 netcdf jsonc sqlite proj gdal metis ; do 
    load_or_install_module $m
done

cmdline="module -t avail netcdf hdf5 metis petsc gdal jsonc proj sqlite"
echo "What is available?"
echo $cmdline
eval $cmdline
