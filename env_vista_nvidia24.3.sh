##
## Basic setup
##
module purge
export MODULEPATH=/opt/apps/modulefiles:/opt/apps/cuda12_4/modulefiles
## :/opt/apps/nvidia24/openmpi_5/modulefiles:/opt/apps/nvidia24/modulefiles

##
## for openmpi
##
export MODULEPATH=/opt/apps/nvidia24/modulefiles

##
## My modules
##
export MODULEPATH=${WORK}/modulefiles/Core:${MODULEPATH}
module load nvhpc-hpcx/24.3
