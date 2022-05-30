% open file and convert to double format
imagefl = imread('Y:\cbarnett\year_greenland_vv\greenland_vel_mosaic_2009_2010_2012_2013.tif');
imagefl_double = double(imagefl);

% obtain spatial reference info
worldfile = getworldfilename('Y:\cbarnett\year_greenland_vv\greenland_vel_mosaic_2009_2010_2012_2013.tif');
SRO = worldfileread(worldfile,'geographic',size(imagefl_double));
% filename and extension for exporting
fn = 'greenland_vel_mosaic500_2009_2010_2012_2013';
fn_ext = '.tif';
full_fn = [fn fn_ext];
% export
geotiffwrite(full_fn, imagefl_double, SRO);
