#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Specify number of cores
#$ -pe omp 8
# Specify memory per core
#$ -l mem_per_core=8G

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=240:00:00

# Name of job
#$ -N calc_mets

# Combine output/error files into single file
#$ -j y

module load matlab/2022b
matlab -nodisplay -r "calc_metrics_main; exit"
