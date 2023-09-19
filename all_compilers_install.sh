#!/bin/bash

if [ $# -lt 1 ] ; then 
    echo "Usage: $0 package" && exit 1
fi

package=$1
for compiler in intel19 intel23 gcc9 gcc13 ; do
    ( \
	cd $package \
	&& make default_install \
    )
done
