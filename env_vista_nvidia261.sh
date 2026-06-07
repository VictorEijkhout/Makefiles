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
NV_VERSION=26.1
module unload nvidia openmpi cuda 2>/dev/null
export MODULEPATH=${WORK}/modulefiles/Core\
:${WORK}/modulefiles/Compiler/nvidia/${NV_VERSION}\
:${MODULEPATH}
module -t load nvidia/${NV_VERSION}
module -t list 2>&1 \
    | awk '{v=v " " $1} END { print "Modules: " v }'
#
echo -e "\nPossible openmpi:"
module -t avail openmpi
OMPI_VERSION=5.0.9
module -t load openmpi/${OMPI_VERSION}
module -t list 2>&1 \
    | awk '{v=v " " $1} END { print "Modules: " v }'
#
echo -n "\nPossible cuda:"
module -t avail cuda
cudaversion=13.1
module -t load cuda/${cudaversion}
MODULEPATH=\
${WORK}/modulefiles/Compiler/nvidia/${NV_VERSION}\
:${MODULEPATH}
echo "now available openmpi:"
module -t avail openmpi
module -t load openmpi/${OMPI_VERSION}

echo "Load nvpl/math"
module -t load nvpl nvidia_math
# #module load nvhpc-hpcx
# module -t list 2>&1 \
#     | awk '{v=v " " $1} END { print "Modules: " v }'

# ##
# ## python hack
# ##
# export TACC_PYTHON_DIR=/opt/apps/gcc14/cuda12/python3/3.11.8
# export TACC_PYTHON_BIN=/opt/apps/gcc14/cuda12/python3/3.11.8/bin
# export PATH=${PATH}:${TACC_PYTHON_BIN}
# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${TACC_PYTHON_DIR}/lib
echo -e "Done setting up nvidia 26 environment\n"
