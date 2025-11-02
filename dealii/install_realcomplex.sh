#!/bin/bash

function usage () {
    echo "Usage: $0 [ -h ] [ -j 123 ] [ -c : only complex ] [ -r : only real ]"
}

jcount=8
real=1
complex=1

while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 0
    elif [ $1 = "-c" ] ; then
	complex=1 && real=0 && echo " .. only installing complex" && shift
    elif [ $1 = "-r" ] ; then
	real=1 && complex=0 && echo " .. only installing real" && shift
    elif [ $1 = "-j" ] ; then
	shift && jcount=$1 && shift
	echo " .. using jcount=$jcount"
    else
	echo "Unknown option <<$1>>" && exit 1
    fi
done

module load $( make listmodules | grep -v petsc | grep -v mkl )

if [ $real -gt 0 ] ; then
    module load petsc/3.24
    make real JCOUNT=$jcount
fi
if [ $complex -gt 0 ] ; then
    module load petsc/3.24-complex
    make complex JCOUNT=$jcount
fi



