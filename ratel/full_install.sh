#!/bin/bash

cd ~/Software/petsc
make pull 
make biginstall PACKAGEVERSION=git CUSTOMEXT=ratel JCOUNT=8

cd ../libceed
make pull
make default_install JCOUNT=8

cd ../ratel
make pull 
make default_install JCOUNT=8
