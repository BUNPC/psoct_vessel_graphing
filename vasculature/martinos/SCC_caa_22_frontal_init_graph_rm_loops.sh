#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Request a whole node with 28 cores and at least 384 GB of RAM.
# Specify number of cores
#$ -pe omp 8
# Specify memory per core
#$ -l mem_per_core=8G

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=240:00:00

# Name of job
#$ -N caa_22f_rmloops

# Combine output/error files into single file
#$ -j y

module load matlab/2022b
matlab -nodisplay -r "caa_22_frontal_remove_loops_main; exit"