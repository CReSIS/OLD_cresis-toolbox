%% Script geotiff_io.m

% Takes input georeference TIFF image file and does the following:
%      1. Put the pixel values into a matrix, find the geographic reference
%         values and the bounding box around the GeoTIFF.
%      2. Find the projection information of the same GeoTIFF.
%      3. Write user-defined calculations.
%      4. Convert the pixel values to unsigned 16-bit intergers (necessary
%         for geotiffwrite.)
%      5. Write modified pixel vales to new GeoTIFF. 
%
% Known issues:
%      1. geotiffwrite only exists in MATLAB 2011b or newer. As of June
%         2012, this does not work in the Windows version of MATLAB.
%
% Date: 21 June 2012
% Author: Steven Foga
%
% See also: geotiffread.m, geotiffwrite.m, geotiffinfo.m, im2uint16.m,
%           mapshow.m

%% Part 1: User Inputs

% Desired Output File
filename = '/cresis/scratch1/cbranecky/GIS/Hq_11_calc.tif';
% filename = 'Z:\sfoga\Hq_11_calc.tif');

% Read in GeoTIFF
[A REFMAT BBOX] = geotiffread('/cresis/scratch1/cbranecky/GIS/Hq_11.tif');
% [A REFMAT BBOX] = geotiffread('Z:\sfoga\Hq_11.tif');

% Get GeoTIFF Projection Information
info = geotiffinfo('/cresis/scratch1/cbranecky/GIS/Hq_11.tif');
% info = geotiffinfo('Z:\sfoga\Hq_11.tif');


%% Part 2: <Your Calculations Here>

% Example: Root transform on GeoTIFF
% A = sqrt(A);

%% Part 3: Automated Output
% Convert raster values from single to uint16
A = im2uint16(A);

geotiffwrite(filename,A,REFMAT,'GeoKeyDirectoryTag', ...
    info.GeoTIFFTags.GeoKeyDirectoryTag);
