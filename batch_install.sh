#!/bin/bash

package=$1

if [ $# -eq 2 ] ; then
    queue=$2
else
    queue=normal
fi

cat >install.slurm <<EOF
#!/bin/bash

#SBATCH -J install${package}           # Job name
#SBATCH -o install${package}.o%j       # Name of stdout output file
#SBATCH -e install${package}.o%j       # Name of stderr error file
#SBATCH -p ${queue}             # Queue (partition) name
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 4:0:0
##SBATCH --mail-user=myname@myschool.edu
##SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A A-ccsc       # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...

cd ${package}
make default_install JCOUNT=20
EOF

sbatch install.slurm
