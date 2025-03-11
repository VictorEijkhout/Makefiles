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
