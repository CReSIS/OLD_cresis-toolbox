#!/bin/bash
# bash script for creating the image browser for release
# First argument: Location of the cresis-toolbox directory
# Second argument: Location of where the image_browser directory will be created
#
# Author: John Paden

if (( `echo $1 | wc -c` < 2 ))
then
  echo "First argument is empty, should be cresis-toolbox location"
  exit
fi
if (( `echo $2 | wc -c` < 2 ))
then
  echo "Second argument is empty, should be output location"
  exit
fi

echo Cresis toolbox at: $1
echo Creating image browser at: $2/image_browser

rm -rf $2/image_browser
mkdir $2/image_browser
mkdir $2/image_browser/display
mkdir $2/image_browser/gps
mkdir $2/image_browser/gui
mkdir $2/image_browser/picker
mkdir $2/image_browser/processing
mkdir $2/image_browser/utility

cp $1/display/lp.m $2/image_browser/display/
cp $1/display/local_detrend.m $2/image_browser/display/
cp $1/gps/epoch_to_datenum.m $2/image_browser/gps/
cp $1/gui/get_mouse_info.m $2/image_browser/gui/
cp -R $1/gui/table $2/image_browser/gui/
cp $1/picker/picker*.m $2/image_browser/picker/
cp $1/picker/run_picker_release.m $2/image_browser/run_picker_release.m
cp $1/processing/track_layer.m $2/image_browser/processing/
cp $1/utility/get_filenames.m $2/image_browser/utility/
cp $1/ct_support/ct_filename_tmp.m $2/image_browser/utility/
cp $1/utility/merge_structs.m $2/image_browser/utility/
cp $1/utility/physical_constants.m $2/image_browser/utility/
cp $1/utility/frame_id_comp.m $2/image_browser/utility/


