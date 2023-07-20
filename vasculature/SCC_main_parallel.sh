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
#$ -l h_rt=48:00:00

# Name of job
#$ -N segment

# Combine output/error files into single file
#$ -j y

# Declare array job (create new job for each)
#$ -t 1-2

echo "Starting task number $SGE_TASK_ID"
module load matlab/2022b
matlab -nodisplay -batch "subid_idx='$SGE_TASK_ID'; psoct_vessel_segmentation_main"


#for subject_id in $(seq $SGE_TASK_ID $(( $SGE_TASK_ID + $SGE_TASK_STEPSIZE - 1 )) ); do
#    echo "Subtask $subtask_id of task $SGE_TASK_ID"
#    ./myprog $subject_id
#done
