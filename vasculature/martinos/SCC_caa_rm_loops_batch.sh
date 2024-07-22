#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Request cores and memory
# Specify number of cores
#$ -pe omp 8
# Specify memory per core
#$ -l mem_per_core=8G

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=240:00:00

# Name of job
#$ -N caa_rmloop

# Combine output/error files into single file
#$ -j y

# Batch array
#$ -t 1-4

echo "Starting task number $SGE_TASK_ID"
module load matlab/2022b
matlab -nodisplay -batch caa_remove_loops_main $SGE_TASK_ID