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
module load nvpl
