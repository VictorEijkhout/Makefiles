##
## Basic setup
##
module purge
module reset

##
## My modules
##
export MODULEPATH=${WORK}/modulefiles/Core:${MODULEPATH}
module load nvidia/24.9
module load openmpi/5.0.5_nvc249
module load nvpl nvidia_math
#module load nvhpc-hpcx

##
## python hack
##
export TACC_PYTHON_DIR=/opt/apps/gcc14/cuda12/python3/3.11.8
export PATH=${PATH}:/opt/apps/gcc14/cuda12/python3/3.11.8/bin/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${TACC_PYTHON_DIR}/lib
