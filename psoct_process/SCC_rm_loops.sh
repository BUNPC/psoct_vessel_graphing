#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Request a whole node with 16 cores and at least 256 GB of RAM.
# Specify number of cores
#$ -pe omp 16
# Specify memory per core
#$ -l mem_per_core=16G

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=24:00:00

# Name of job
#$ -N segment

# Combine output/error files into single file
#$ -j y

module load matlab/2022b
matlab -nodisplay -r "rm_loops_main; exit"

