#!/bin/bash

# Add some debugging information to stdout
echo "Running cresis-toolbox cluster_job.{sh,m}"
id
hostname
declare
date
umask 000
pwd

# Turn on MCR debugging information (goes to stderr)
export  MCR_CACHE_VERBOSE=1

# This is a Matlab compiled .m file which calls the specific cresis-toolbox command
# specified in the environment variables (see worker_task.m for details)

echo "Calling run_cluster_job.sh"
$MATLAB_CLUSTER_PATH/run_cluster_job.sh $MATLAB_MCR_PATH
echo "Returned from run_cluster_job.sh"

date

