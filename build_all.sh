#!/bin/bash
################################################################
####
#### Build all package, respecting module dependencies
#### - this assumes that `all_packages/all_ladder.sh' exists
#### - which is created by running `build_ladder.sh'
####
################################################################

setx=0
jcount=4
list=
packages=0

phdf5_dir=hdf5
phdf5_tgt=par
PACKAGEOPTIONS_hdf5="HDFFORTRAN=OFF"

petsc_commandline="./install_big.sh -f"

# these are module names
cd all_packages
ladderfile="all_ladder.sh"
if [ ! -f "${ladderfile}" ] ; then
    echo "ERROR could not find ladder file <<${ladderfile}>>"
    echo " first run <<build_ladder.sh>>"
    exit 1
fi
../ladder.sh $* "${ladderfile}"


