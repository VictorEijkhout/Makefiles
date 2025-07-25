## where do the spaces in this come from?
numladder="\
    $( i=1 && for l in ${ladder} ; do echo $i,$l && i=$(( i+1 )) ; done ) \
    "
export ladderlength=$( echo "${numladder}" | wc -w )

function usage() {
    echo "Usage: $0 [ -h ] [ -x : setx ] [ -l : list ] [ -t : trace ]"
    echo "    [ -j nnn (default: ${jcount}) ] "
    echo "    nnn / all"
    echo "where nnn:"
    for ixy in ${numladder} ; do
	n=${ixy%%,*}
	xy=${ixy#*,}
	x=${xy%,*}
	y=${xy#*,}
	echo " $n : $x / $y"
    done
    echo " all : 1--${n}"
}
if [ $# -eq 0 ] ; then 
    usage
    exit 0
fi

function parse_numpacver {
    numpacver="$1"
    # number in the list
    num=${numpacver%%,*}
    # package and optional version
    pacver=${numpacver#*,}
    # package
    package=${pacver%,*}
    # directory, could be package name or different
    eval packagedir=\${${package}_dir}
    if [ -z "${packagedir}" ] ; then packagedir=${package} ; fi
    # version could be empty
    loadversion=$( pacver="${pacver}," && echo ${pacver#*,} | tr -d ',' )
    if [ ! -z $loadversion ] ; then 
	loadversion=$( cd ../${packagedir} && make --no-print-directory loadversion )
    fi
}

function module_avail {
    # $1=package $2=test version $3=install version
    if [ -z "$2" ] ; then 
	testversion="$1/$3"
    else
	testversion="$1/$2"
    fi
    avail=$( module -t avail ${testversion} 2>&1 )
    if [ -z "${avail}" ] ; then 
	echo "MISSING: $testversion:"
    else
	echo "Available: $testversion"
	echo "${avail}"
    fi
}

function module_install {
    package="$1" fullversion="$2"
    echo " .. installing" && echo
    eval packagecmdline=\${${package}_commandline}
    if [ ! -z "${packagecmdline}" ] ; then
	pushd ../${packagedir} 
	eval ${packagecmdline} 
	popd
    else
	eval packagetgt=\${${package}_tgt}
	if [ -z "${packagetgt}" ] ; then packagetgt=default_install ; fi 
	pushd ../${packagedir} \
	    && make ${packagetgt} public \
		    $( if [ "${trace}" = "1" ] ; then echo ECHO=1 ; fi ) \
		    JCOUNT=${jcount} versionspec="PACKAGEVERSION=$fullversion"
	popd
    fi
}

##
## Commandline options
##
while [ $# -gt 0 ] ; do
    if [ "$1" = "-h" ] ; then 
	usage
	exit 0
    elif [ "$1" = "-l" ] ; then 
	list=1; shift
    elif [ "$1" = "-x" ] ; then 
	setx=1; shift
    elif [ "$1" = "-j" ] ; then 
	shift; jcount=$1; shift
    elif [ "$1" = "-t" ] ; then 
	trace=1; shift
    else
	break
    fi
done

if [ $# -eq 0 ] ; then
    # this can only happen in list mode
    if [ -z "${list}" ] ; then
	usage && exit 0
    fi
    packages="$( seq 1 $ladderlength )"
    echo "---------------- listing packages"
elif [ "$1" = "all" ] ; then
    packages="$( seq 1 $ladderlength )"
    echo "---------------- installing packages: ${packages}"
else
    packages=$( echo $* | tr "," " " )
    echo "---------------- installing packages: ${packages}"
fi

if [ $setx -gt 0 ] ; then 
    set -x
fi

##
## Exceptions
##
vista_exclude="fftw2 json sfml"
eval packages_to_exclude=\${${TACC_SYSTEM}_exclude}

##
## Modules without directory
##
phdf5_dir=hdf5
phdf5_tgt=par
parallelnetcdf_dir=netcdf
parallelnetcdf_tgt=par

##
## Log
##
ladderlog=ladder_${TACC_FAMILY_COMPILER}${TACC_FAMILY_COMPILER_VERSION}.log
( \
  echo "================ Starting installation with modules:" \
      && module -t list 2>&1 | sort | tr '\n' ' ' && echo \
  ) | tee ${ladderlog}

##
## Install loop
##
for m in \
    $( for p in ${packages} ; do 
        if [[ $p = *\-* ]] ; then 
	    # expand rangers
	    echo $( seq $( echo $p | cut -d '-' -f 1 ) $( echo $p | cut -d '-' -f 2 ) )
	elif [[ $p != \~* ]] ; then
	    # list number if not negated
	    echo $p
	fi
       done ) ; do
    echo "================================================================"
    echo "==== Load or build: $m"
    echo "================================================================"
    for numpacver in ${numladder} \
	       ; do \
	parse_numpacver "${numpacver}"
	echo "================"
	echo "Package $num: $package version $loadversion"
	if [ ! -z "${list}" ] ; then 
	    module_avail "$package" "$loadversion" "$fullversion"
	else
	    # system-dependent exclude
	    # otherwise install or load
	    if [ $m -eq $num ] ; then 
		if [[ ${packages_to_exclude} != *${package}* ]] ; then 
		    module_install "${package}" "${fullversion}"
		fi
		# go to next installable
		break
	    else
		echo "Loading package/version = $package/$loadversion"
		if [[ ${packages_to_exclude} = *${package}* ]] ; then
		    echo " .. ignore load" && continue
		fi
		module load $package/$loadversion
		if [ $? -ne 0 ] ; then echo "Could not load $package" && exit 1 ; fi
		module -t list $package/$loadversion 
		PACKAGE=$( echo ${package} | tr a-z A-Z )
		eval packagedir=\${TACC_${PACKAGE}_DIR}
		if [ -z "${packagedir}" ] ; then 
		    echo "ERROR null packagedir variable for package=${package}"
		    exit 1
		elif [ ! -d "${packagedir}" ] ; then 
		    echo "ERROR no such packagedir: ${packagedir}"
		fi
		eval echo " .. loaded ${package}/${loadversion} at ${packagedir}"
	    fi
	fi
    done 
    if [ ! -z "${list}" ] ; then break ; fi
done 2>&1 | tee -a ${ladderlog}

echo && echo "See: ${ladderlog}" && echo


