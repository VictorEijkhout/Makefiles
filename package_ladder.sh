#!/bin/bash

function usage () {
    echo "Usage: $0 [ -h ] [ -j 123 ] ladderfile"
}

if [ $# -eq 0 ] ; then usage && exit 0 ; fi

while [ $# -gt 1 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 1
    elif [ $1 = "-j" ] ; then
	shift && jcount="JCOUNT=$1" && shift
    else
	echo "Unknown option <<$1>>" && exit 1
    fi
done
ladderfile=$1
if [ ! -f ${ladderfile} ] ; then
    echo "ERROR can not find ladder file <<$ladderfile>>"
    exit 1
fi

phdf5_dir=hdf5
phdf5_tgt=par

logfile=$(pwd)/ladder_${TACC_FAMILY_COMPILER}_${TACC_FAMILY_COMPILER_VERSION}.log
for package in $( cat ladder.txt ) ; do
    module -t load ${package} 2>/dev/null
    if [ $? -eq 0 ] ; then
	echo
	echo "Package <<${package}>> successfully loaded:"
	module -t list ${package}
	echo 
    else
	if [ -d ../${package} ] ; then
	    echo "${package}: making default target"
	    pushd ../${package}
	    make default_install ${jcount}
	    popd
	else
	    eval dir=\${${package}_dir}
	    eval tgt=\${${package}_tgt}
	    if [ -z "${dir}" ] ; then
		echo "ERROR no directory given for <<${package}>>" && exit 1
	    fi
	    if [ -z "${tgt}" ] ; then
		echo "ERROR no target given for <<${package}>>" && exit 1
	    fi
	    pushd ../${dir}
	    make ${tgt} ${jcount}
	    popd
	fi
	module load ${package}
    fi
done 2>&1 | tee ${logfile}

echo
echo "See: ${logfile}"
echo
