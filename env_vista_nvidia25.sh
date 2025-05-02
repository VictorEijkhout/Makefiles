##
## Basic setup
##
module purge
module reset

##
## My modules
##
export MODULEPATH=${WORK}/modulefiles/Core:${MODULEPATH}
module load nvidia/25.3
module load openmpi/5.0.7
module load nvpl nvidia_math
#module load nvhpc-hpcx

##
## python hack
##
export TACC_PYTHON_DIR=/opt/apps/gcc14/cuda12/python3/3.11.8
export TACC_PYTHON_BIN=/opt/apps/gcc14/cuda12/python3/3.11.8/bin
export PATH=${PATH}:${TACC_PYTHON_BIN}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${TACC_PYTHON_DIR}/lib
