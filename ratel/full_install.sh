#!/bin/bash
################################################################
#### install all software, given externally loaded compiler
################################################################

cd ~/Software/petsc
make pull 
make biginstall PACKAGEVERSION=git CUSTOMEXT=ratel JCOUNT=8
if [ $? -gt 0 ] ; then exit 1; fi
module load petsc/git-ratel

cd ../libceed
make pull
make default_install JCOUNT=8
if [ $? -gt 0 ] ; then exit 1; fi
module load libceed/git

cd ../ratel
make pull 
make default_install JCOUNT=8
if [ $? -gt 0 ] ; then exit 1; fi
module load ratel

cd ../tau
make configure build JCOUNT=8
if [ $? -gt 0 ] ; then exit 1; fi
module load tau
