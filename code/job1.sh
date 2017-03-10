#!/bin/bash\n
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --job-name=flintstones
mpirun ./flintstone
