#!/bin/bash

function usage () {
    echo "Usage: $1 [ -h ] [ --install big/small]"
}

export INSTALLATION=big
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit
    elif [ $1 == "--install" ] ; then
	shift && INSTALLATION=$1 && shift
    else
	echo "Unknown option <<$1>>" && usage && exit 1
    fi
done

for i in 32 46 ; do
    export INTSIZE=$i
    for s in real complex ; do
	export SCALAR=$s
	for d in OFF ON ; do
	    export DEBUG=$d
	    SCRIPTSDIR=${PWD}/mpm_scripts_$s$i$d \
		      mpm.py install
	done
    done
done
