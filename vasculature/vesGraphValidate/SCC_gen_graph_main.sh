#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Request a whole node with 28 cores and at least 384 GB of RAM.
# Specify number of cores
#$ -pe omp 16
# Specify memory per core
#$ -l mem_per_core=16G

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=24:00:00

# Name of job
#$ -N gen_graph

# Combine output/error files into single file
#$ -j y

# Declare array job (create new job for each)
# Iterate over all subjects and use smallest sigma (1-51:3)
# Iterate over all subjects and use all sigmas (1-51)
#$ -t 1-15

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="

module load matlab/2022b
matlab -nodisplay -singleCompThread -r "gen_graph_main; exit"

