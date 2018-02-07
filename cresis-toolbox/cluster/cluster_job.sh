#!/bin/bash

# =========================================================================
# SETUP
# =========================================================================

# Add some debugging information to stdout
echo "cluster_job.{sh,m} Start" `whoami` `hostname` "("`date`")"
#declare
#pwd

# Make sure file permissions will be set so everyone can read
umask 000

# Turn on MCR debugging information (goes to stderr)
export  MCR_CACHE_VERBOSE=1

# Get this scripts (the parent) process info
parent_pid=$$
#echo Parent: $parent_pid

# =========================================================================
# ATTEMPT TO RUN CLUSTER_JOB.M UP TO MAX_ATTEMPTS TIMES
# =========================================================================
max_attempts=3
attempt=1

while (( attempt < max_attempts ))
do

  # Start child process
  # Run run_cluster_job.sh (runs Matlab compiled cluster_job.m)
  $MATLAB_CLUSTER_PATH/run_cluster_job.sh $MATLAB_MCR_PATH &

  child_pid=$!
  #echo Child: $child_pid

  # Wait for child process to finish
  child_proc=`ps -eo ppid,pid | sed -n "/^\s*$parent_pid\s*$child_pid/p" | awk '{print $2}'`
  max_mem=0
  max_cpu=0
  while [[ ! -z $child_proc ]]
  do
    mem=`ps -eo ppid,pid,rss | sed -n "/^\s*$parent_pid\s*$child_pid/p" | awk '{print $3}'`
    #echo Mem: $mem
    if (( mem > max_mem ))
    then
      max_mem=$mem
    fi
    cpu=`ps -eo ppid,pid,cputime | sed -n "/^\s*$parent_pid\s*$child_pid/p" | awk '{print $3}'`
    #echo CPU: $cpu
    if [[ ! -z $cpu ]]
    then
      max_cpu=$cpu
    fi

    # Update maximum memory and maximum CPU

    child_proc=`ps -eo ppid,pid | sed -n "/^\s*$parent_pid\s*$child_pid/p" | awk '{print $2}'`
    sleep 1
  done
  echo Max Mem: $max_mem
  echo Max CPU: $max_cpu

  # Get return value from child process
  wait $child_pid
  return_value=$?
  echo Return value: $return_value

  if [ $return_value -ne 255 ]
  then
    break
  fi
  ((attempt++))

done

echo "  cluster_job.{sh,m} Done" `whoami` `hostname` "("`date`")"
date

exit $return_value
