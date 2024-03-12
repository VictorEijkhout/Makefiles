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
pversion=3.20.4
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
archs=archs-${pversion}
rm -f $archs && touch $archs
for debug in 0 1 ; do 
    for int in 32 64 ; do 
	for precision in single double ; do
	    for scalar in real complex ; do 
		export DEBUG=${debug}
		export INT=${int}
		export PRECISION=${precision}
		export SCALAR=${scalar} 
		arch=$( make --no-print-directory petscshortarch )
		echo && echo "Installing big for arch=${arch}" && echo
		if [ -z "${arch}" ] ; then 
		    echo vanilla >> ${archs}
		else
		    echo ${arch} >> ${archs}
		fi
		./install_big.sh -j ${jcount} \
		    $( if [ ${p4p} -eq 0 ] ; then echo "-4" ; fi ) \
		    $( if [ ${cuda} -eq 1 ] ; then echo "-c" ; fi )
	    done
	done
    done
done

export NOHDF5=1
for debug in 0 1 ; do 
    export INT=32
    export PRECISION=double
    export SCALAR=real 
    export DEBUG=${debug}
    arch=$( make --no-print-directory petscshortarch )
    echo && echo "Installing big for arch=${arch}" && echo
    echo ${arch} >> ${archs}
    ./install_big.sh -j ${jcount} \
    	-5 \
    	$( if [ ${p4p} -eq 0 ] ; then echo "-4" ; fi ) \
    	$( if [ ${cuda} -eq 1 ] ; then echo "-c" ; fi )
done
export NOHDF5=0

export FORTRAN=08
for debug in 0 1 ; do 
    for scalar in real complex ; do 
	export INT=32
	export PRECISION=double
	export SCALAR=${scalar}
	export DEBUG=${debug}
	arch=$( make --no-print-directory petscshortarch )
	echo && echo "Installing big for arch=${arch}" && echo
	echo ${arch} >> ${archs}
	./install_big.sh -j ${jcount} \
	    -8 \
	    $( if [ ${p4p} -eq 0 ] ; then echo "-4" ; fi ) \
	    $( if [ ${cuda} -eq 1 ] ; then echo "-c" ; fi )
    done
done

echo && echo "done archs: $( cat ${archs} | tr '\n' ' ' )" && echo
