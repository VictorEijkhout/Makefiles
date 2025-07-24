#!/bin/bash

setx=0
jcount=4
list=
packages=0

phdf5_dir=hdf5
phdf5_tgt=par
PACKAGEOPTIONS_hdf5="HDFFORTRAN=OFF"

petsc_commandline="./install_big.sh -f"

# these are module names
ladder="$( cat all_ladder.sh )"

source ../ladder.sh

