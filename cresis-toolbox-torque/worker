#!/bin/bash

# Add some debugging information to stdout
echo "Running cresis-toolbox Worker"
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

echo "Calling run_worker_task.sh"
$MATLAB_TORQUE_PATH/run_worker_task.sh $MATLAB_MCR_PATH
echo "Returned from run_worker_task.sh"

date

