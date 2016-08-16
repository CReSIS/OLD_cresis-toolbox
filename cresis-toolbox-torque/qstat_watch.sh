#!/bin/bash

user=`whoami`
TMP_FILE=/tmp/qstat_watch.$USER.tmp
LOG_FILE=/tmp/qstat_watch.$USER.log

printf "Writing to file %s\n" $LOG_FILE

while (( 1 ))
do
  datestr=`date`
  printf "=============== $datestr ================\n"
  printf "=============== $datestr ================\n" >>$LOG_FILE
  #TMP_FILE=/tmp/qstat_watch.$$.tmp
  qstat -u $USER >$TMP_FILE
  queued=`grep " Q " $TMP_FILE | wc -l`
  running=`grep " R " $TMP_FILE | wc -l`
  printf "There are %d queued and %d running\n" $queued $running
  printf "There are %d queued and %d running\n" $queued $running >>$LOG_FILE
  sed -n -e 3,4p $TMP_FILE
  sed -n -e 3,4p $TMP_FILE >>$LOG_FILE
  grep " Q " $TMP_FILE | head -1
  grep " Q " $TMP_FILE | head -1 >>$LOG_FILE
  grep " R " $TMP_FILE | head -1
  grep " R " $TMP_FILE | head -1 >>$LOG_FILE
  #awk '{ print $11}' /tmp/qstat_watch.tmp
  sleep 20
done

