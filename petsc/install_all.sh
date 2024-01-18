#!/bin/bash

function usage() {
    echo "Usage: $0 [ -h ] "
    echo "    [ -v (default: ${pversion} ]"
    echo "    [ -j jpar (default: ${jcount}) ]"
    echo "    [ -c : cuda build ]"
    echo "    [ -4 : skip petsc/slepc4py ]"
}

cuda=0
jcount=6
p4p=1
pversion=3.20.3
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 0
    elif [ $1 = "-4" ] ; then 
	p4p=0 && shift
    elif [ $1 = "-c" ] ; then 
	cuda=1 && shift
    elif [ $1 = "-j" ] ; then
	shift && jcount=$1 && shift
    elif [ $1  = "-v" ] ; then 
	shift && pversion=$1 && shift
    else 
	echo "Error: unknown option <<$1>>" && exit 1
    fi
done

set -e
for debug in 0 1 ; do 
    for int in 32 64 ; do 
	for precision in single double ; do
	    for scalar in real complex ; do 
		export DEBUG=${debug}
		export INT=${int}
		export PRECISION=${precision}
		export SCALAR=${scalar} 
		./install_big.sh -j ${jcount} \
		    $( if [ ${p4p} -eq 0 ] ; then echo "-4" ; fi ) \
		    $( if [ ${cuda} -eq 1 ] ; then echo "-c" ; fi )
	    done
	done
    done
done

