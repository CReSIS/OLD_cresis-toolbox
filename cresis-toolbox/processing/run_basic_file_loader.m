% script run_basic_file_loader
%
% Example script for running basic_file_loader.m
%
% Author: John Paden

[param,defaults] = default_radar_params_2016_Greenland_Polar6_mcords;

[data,fn,settings,default,gps,hdr,pc_param] = basic_file_loader(param,defaults);
