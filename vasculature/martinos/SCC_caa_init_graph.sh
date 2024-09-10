#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Specify number of cores
#$ -pe omp 16
# Specify memory per core
#$ -l mem_per_core=16G

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=240:00:00

# Name of job
#$ -N caa_rmloop

# Combine output/error files into single file
#$ -j y

# Batch array
#$ -t 3

echo "Starting task number $SGE_TASK_ID"
module load matlab/2022b
matlab -nodisplay -batch caa_init_graph $SGE_TASK_ID