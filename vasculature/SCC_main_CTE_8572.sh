#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Request a whole node with 16 cores and at least 128 GB of RAM.
# Specify number of cores
#$ -pe omp 16

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=24:00:00

# Name of job
#$ -N segment_CTE_8572

# Combine output/error files into single file
#$ -j y

module load matlab/2022b
matlab -nodisplay -singleCompThread -r "psoct_vessel_segmentation_main_AD_10382; exit"

