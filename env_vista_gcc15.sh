## module load lmod
module -t reset 2>/dev/null
module -t use ${WORK}/modulefiles/Core
module -t load gcc/15.1 2>/dev/null
module -t load openmpi/5.0.8 2>/dev/null
module -t load cuda nvpl nvidia_math 2>/dev/null
module -t list 2>&1 | sort | tr '\n' ' '
echo



