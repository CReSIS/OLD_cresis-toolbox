#!/bin/bash

echo $1
#TAG=1chop
TAG=2chop_picks
COUNT=`echo images/${1}/*_${TAG}* ${1}_combined.png | wc -w`
echo $COUNT
echo images/${1}/*_${TAG}* ${1}_picks_combined.png |  xargs montage +frame +shadow +label -geometry +0+0 -tile ${COUNT}x1 -depth 8 -type TrueColor


