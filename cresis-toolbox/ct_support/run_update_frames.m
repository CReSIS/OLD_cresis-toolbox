% script run_update_frames
%
% Script for running "update_frames.m"

clear param;
param.radar_name = 'snow';
% param.day_seg = '20090425_04';
% param.day_seg = '20100421_03';
param.day_seg = '20120402_01';

% param.season_name = '2009_Greenland_P3';
% param.image_out_dir = {{'/cresis/snfs1/dataproducts/ct_data/snow/2009_Greenland_P3/CSARP_post/images/'}, ...
%   {'/cresis/snfs1/dataproducts/public/data/temp/for_Snow_Workshop/2009_Greenland_P3/images/'}};
% param.mat_out_dir = {{'/cresis/snfs1/dataproducts/ct_data/snow/2009_Greenland_P3/CSARP_post/CSARP_deconv/'}, ...
%   {'/cresis/snfs1/dataproducts/public/data/temp/for_Snow_Workshop/2009_Greenland_P3/CSARP_deconv/'}};
% param.image_out_dir = {{'/cresis/snfs1/dataproducts/ct_data/snow/2009_Greenland_P3/CSARP_post/images/'}};
% param.mat_out_dir = {{'/cresis/snfs1/dataproducts/ct_data/snow/2009_Greenland_P3/CSARP_post/CSARP_deconv/'}};

% param.season_name = '2010_Greenland_DC8';
% param.image_out_dir = {{'/cresis/snfs1/dataproducts/public/data/snow/2010_Greenland_DC8/images/'}, ...
%   {'/cresis/snfs1/dataproducts/public/data/temp/for_Snow_Workshop/2010_Greenland_DC8/images/'}};
% param.mat_out_dir = {{'/cresis/snfs1/dataproducts/public/data/snow/2010_Greenland_DC8/CSARP_deconv/'}, ...
%   {'/cresis/snfs1/dataproducts/public/data/temp/for_Snow_Workshop/2010_Greenland_DC8/CSARP_deconv/'}};
% param.image_out_dir = {{'/cresis/snfs1/dataproducts/public/data/snow/2010_Greenland_DC8/images/'}};
% param.mat_out_dir = {{'/cresis/snfs1/dataproducts/public/data/snow/2010_Greenland_DC8/CSARP_deconv/'}};

% param.season_name = '2011_Greenland_P3';
% param.image_out_dir = {{'/cresis/snfs1/dataproducts/public/data/snow/2011_Greenland_P3/images/'}, ...
%   {'/cresis/snfs1/dataproducts/public/data/temp/for_Snow_Workshop/2011_Greenland_P3/images/'}};
% param.mat_out_dir = {{'/cresis/snfs1/dataproducts/public/data/snow/2011_Greenland_P3/CSARP_deconv/'}, ...
%   {'/cresis/snfs1/dataproducts/public/data/temp/for_Snow_Workshop/2011_Greenland_P3/CSARP_deconv/'}};
% param.image_out_dir = {{'/cresis/snfs1/dataproducts/public/data/snow/2011_Greenland_P3/images/'}};
% param.mat_out_dir = {{'/cresis/snfs1/dataproducts/public/data/snow/2011_Greenland_P3/CSARP_deconv/'}};

param.season_name = '2012_Greenland_P3';
% param.image_out_dir = {{'/cresis/snfs1/dataproducts/public/data/snow/2012_Greenland_P3/images/'}, ...
%   {'/cresis/snfs1/dataproducts/public/data/temp/for_Snow_Workshop/2012_Greenland_P3/images/'}};
% param.mat_out_dir = {{'/cresis/snfs1/dataproducts/public/data/snow/2012_Greenland_P3/CSARP_deconv/'}, ...
%   {'/cresis/snfs1/dataproducts/public/data/temp/for_Snow_Workshop/2012_Greenland_P3/CSARP_deconv/'}};
param.image_out_dir = {{'/cresis/snfs1/dataproducts/public/data/snow/2012_Greenland_P3/images/'}};
param.mat_out_dir = {{'/cresis/snfs1/dataproducts/public/data/snow/2012_Greenland_P3/CSARP_deconv/'}};

param.img_type = 'echo.jpg';

% update_field: string containing the field in frames file to update
% ('deconv_wf', 'proc_mode', 'quality', etc)
% update_field = 'deconv_wf'; update_field_type = 'double';
% update_field = 'proc_mode'; update_field_type = 'double';
update_field = 'quality'; update_field_type = 'mask';
update_field_mask = {'turn all masks off','coherent noise','deconvolution artifact', ...
  'raised noise floor/vertical stripes','missing data','no good data','low SNR','unclassified','land or iceberg'};

audio_tone_for_nonzero_nonisnan = true;
%audio_tone_check_code = '~isnan(frames.(update_field)(frm)) && frames.(update_field)(frm) ~= 0';
audio_tone_check_code = 'frames.(update_field)(frm) >= 8';

%update_field_match = [1 2 3];
update_field_match = [];

% FMCW image debugging mode (for snow,kuband radar) 
fmcw_img_debug_mode = true;
noise_time_duration = 45e-9;
noise_time_buffer = 250e-9;
img_sidelobe = -35;
noise_threshold_offset_dB = 3.2;

update_frames;

