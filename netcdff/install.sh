#!/bin/bash

echo && echo "Installing sequential" && echo 
module swap phdf5 hdf5 && module swap pnetcdf netcdf && make seq JCOUNT=10 && module show netcdff

echo && echo "Installing parallel" && echo 
module swap hdf5 phdf5 && module swap netcdf pnetcdf && make par JCOUNT=10 && module show pnetcdff
