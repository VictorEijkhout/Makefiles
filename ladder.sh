## where do the spaces in this come from?
numladder="\
    $( i=1 && for l in ${ladder} ; do echo $i,$l && i=$(( i+1 )) ; done ) \
    "
export ladderlength=$( echo "${numladder}" | wc -w )

compiler=${TACC_FAMILY_COMPILER}
version=${TACC_FAMILY_COMPILER_VERSION}
version=${version%%.*}

function usage() {
    echo "Usage: $0 [ -h ] [ -x : setx ] [ -l : list ] "
    echo "    [ -j nnn (default: ${jcount}) ] "
    echo "    [ -c compiler (default ${TACC_FAMILY_COMPILER} ] "
    echo "    [ -v compiler_version (default ${TACC_FAMILY_COMPILER_VERSION} ] "
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
    version=$( pacver="${pacver}," && echo ${pacver#*,} | tr -d ',' )
    if [ ! -z $version ] ; then 
	eval fullversion=\${${package}_full_version}
	if [ -z "$fullversion" ] ; then fullversion=${version} ; fi
    else 
	# get default from makefile
	fullversion=$( cd ../${packagedir} && make --no-print-directory version )
	version=${fullversion}
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
    eval packagetgt=\${${package}_tgt}
    if [ -z "${packagetgt}" ] ; then packagetgt=default_install ; fi 
    pushd ../${packagedir} \
	&& make ${packagetgt} public \
	JCOUNT=${jcount} versionspec="PACKAGEVERSION=$fullversion"
    popd

}

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
    elif [ "$1" = "-c" ] ; then 
	shift; compiler=$1; shift
    elif [ "$1" = "-v" ] ; then 
	shift; version=$1; shift
    else
	if [ "$1" = "all" ] ; then
	    packages="$( seq 1 $ladderlength )"
	else
	    packages=$1
	fi
	shift
    fi
done

export TACC_FAMILY_COMPILER=${compiler}
export TACC_FAMILY_COMPILER_VERSION=${version}
settings=../env_${TACC_SYSTEM}_${TACC_FAMILY_COMPILER}${TACC_FAMILY_COMPILER_VERSION}.sh
if [ ! -f "${settings}" ] ; then 
    echo "Error: no such settings file <<${settings}>>" && exit 1 ; fi
source ${settings} >/dev/null 2>&1

if [ "${packages}" = "0" -a -z "${list}" ] ; then 
  exit 0
fi

if [ $setx -gt 0 ] ; then 
    set -x
fi

ladderlog=ladder_${TACC_FAMILY_COMPILER}${TACC_FAMILY_COMPILER_VERSION}.log
( \
  echo "================ Starting installation with modules:" \
      && module -t list 2>&1 | sort | tr '\n' ' ' && echo \
  ) | tee ${ladderlog}

if [ -z "${list}" ] ; then
    echo "---------------- installing packages: ${packages}"
else
    echo "---------------- listing packages"
fi
for m in $( echo "${packages}" | tr , ' ' ) ; do
    for numpacver in ${numladder} \
	       ; do \
	parse_numpacver "${numpacver}"
	echo "================"
	echo "Package $num: $package version $version"
	if [ ! -z "${list}" ] ; then 
	    module_avail "$package" "$version" "$fullversion"
	elif [ $m -eq $num ] ; then 
	    module_install "${package}" "${fullversion}"
	    # go to next installable
	    break
	fi
	if [ -z "${list}" ] ; then 
	    module load $package/$version
	    if [ $? -ne 0 ] ; then echo "Could not load $package" && exit 1 ; fi
	    PACKAGE=$( echo ${package} | tr a-z -A-Z )
	    module -t list $package/$version 
	    #module show $package/$version 2>&1 | \grep DIR
	    eval echo " .. loaded ${package}/${version} at \${TACC_${package}_DIR}"
	fi
    done 
done 2>&1 | tee -a ${ladderlog}

echo && echo "See: ${ladderlog}" && echo


