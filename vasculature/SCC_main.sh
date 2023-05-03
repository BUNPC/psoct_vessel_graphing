#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Specify job to run on lab's buy-in computers on SCC
#$ -l buyin


# Request a whole node with 28 cores and at least 384 GB of RAM.
# Specify number of cores
#$ -pe omp 16
# Specify memory per core
#$ -l mem_per_core=13G

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=24:00:00

# Name of job
#$ -N ves_seg_graph

# Combine output/error files into single file
#$ -j y

module load matlab/2020b
matlab -nodisplay -singleCompThread -r "psoct_vessel_segmentation_main; exit"

