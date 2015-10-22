#!/bin/bash

strlen1=`echo $1 | wc -c`
strlen2=`echo $2 | wc -c`

if (( strlen1 <= 1 || strlen2 <=1 ))
then
  echo Script for copying just the essential files of cresis-toolbox-torque
  echo Argument 1: Full path of cresis-toolbox-torque
  echo Argument 2: Destination path of where to put the directory
  echo Example:
  echo " ./copy_cresis-toolbox-torque.sh \`pwd\` ../archive/"
  exit
fi

echo Copying $1 to $2

cmd="mkdir $2/cresis-toolbox-torque"
eval $cmd
cmd="cp -i $1/README_torque.txt $2/cresis-toolbox-torque/"
eval $cmd
cmd="cp -i $1/copy_cresis-toolbox-torque.sh $2/cresis-toolbox-torque/"
eval $cmd
cmd="cp -i $1/worker $2/cresis-toolbox-torque/"
eval $cmd
cmd="cp -i $1/worker_task.m $2/cresis-toolbox-torque/"
eval $cmd

