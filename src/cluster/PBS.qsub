#!/bin/bash
#PBS -S /bin/bash
#PBS -q consort
#PBS -l nodes={{nnodes}}:ppn={{ppnodes}}
#PBS -l walltime={{walltime}}
#PBS -N {{name}}
#PBS -o {{logout}}
#PBS -e {{logerr}}

# Change to working directory
cd {{wd}}

# load all modules required for julia
module load gcc/4.8.1
module load python

echo "calling julia now"
julia -p {{numprocs}} {{run}} 
