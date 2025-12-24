##
## Basic setup
##
module purge
module reset
module -t list 2>&1 \
    | awk '{v=v " " $1} END { print "Initial modules: " v }'

##
## My modules
##
module unload nvidia openmpi cuda 2>/dev/null
export MODULEPATH=${WORK}/modulefiles/Core\
:${WORK}/modulefiles/Compiler/nvidia/25.9\
:${MODULEPATH}
module load nvidia/25.9
module load openmpi/5.0.8
module load cuda/12.9
module load nvpl nvidia_math
#module load nvhpc-hpcx
module -t list 2>&1 \
    | awk '{v=v " " $1} END { print "Modules: " v }'

##
## python hack
##
export TACC_PYTHON_DIR=/opt/apps/gcc14/cuda12/python3/3.11.8
export TACC_PYTHON_BIN=/opt/apps/gcc14/cuda12/python3/3.11.8/bin
export PATH=${PATH}:${TACC_PYTHON_BIN}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${TACC_PYTHON_DIR}/lib
