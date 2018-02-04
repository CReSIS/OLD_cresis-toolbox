#!/bin/bash

# Add some debugging information to stdout
echo "cluster_job.{sh,m} Start" `whoami` `hostname` "("`date`")"
#declare
#pwd

# Make sure file permissions will be set so everyone can read
umask 000

# Turn on MCR debugging information (goes to stderr)
export  MCR_CACHE_VERBOSE=1

# Run run_cluster_job.sh (runs Matlab compiled cluster_job.m)
$MATLAB_CLUSTER_PATH/run_cluster_job.sh $MATLAB_MCR_PATH

echo "  cluster_job.{sh,m} Done" `whoami` `hostname` "("`date`")"
date
