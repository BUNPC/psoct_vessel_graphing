#!/bin/bash -l

# Set SCC project
#$ -P npbssmic

# Send email upon completion
#$ -m ea

# Time limit for job
#$ -l h_rt=1:00:00

# Name of job
#$ -N test_array

# Combine output/error files into single file
#$ -j y

### Declare array job (create new job for each)
# There are 17 subjects and 3 gaussian sigma arrays
# (small, medium, large) for a total of 51 jobs. 
#$ -t 1-2

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="

module load matlab/2022b

# Retrieve contents of line number $SGE_TASK_ID from the text file subject_sigma_index_list
params=`sed -n "${SGE_TASK_ID} p" subject_sigma_index_list.txt`
paramsArray=($params)

# Assign first element of line to subject ID index
subid_idx=${paramsArray[0]}

# Assign second element of line to gaussian sigma array index
gauss_idx=${paramsArray[1]}

matlab -nodisplay -batch "scc_parallel_test" $subid_idx $gauss_idx

