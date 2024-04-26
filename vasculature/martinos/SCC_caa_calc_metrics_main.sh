#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Request a partial node with 4 cores and 64 GB of RAM.
# Specify number of cores
#$ -pe omp 4
# Specify memory per core
#$ -l mem_per_core=16G

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=11:30:00

# Name of job
#$ -N caa_metrics

# Combine output/error files into single file
#$ -j y

module load matlab/2022b
matlab -nodisplay -nodesktop -r "caa_init_graph; exit"
