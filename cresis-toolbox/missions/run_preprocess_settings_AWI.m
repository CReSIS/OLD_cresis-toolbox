% script run_preprocess_AWI.m
%
% Support script for run_preprocess.m

param.config.default = [];

%% MCORDS5 AWI UWB - SINGLE DAY
cur_idx = length(param.config.default)+1;
param.config.default{cur_idx} = default_radar_params_2021_Greenland_Polar5_mcords();
param.config.base_dir{cur_idx} = 'E:\';
param.config.config_folder_names{cur_idx} = '2107300401';
param.config.board_folder_names{cur_idx} = fullfile('2107300401','%b');
param.config.date_strs{cur_idx} = '20210730';

%% SNOW8 SINGLE DAY
% cur_idx = length(param.config.default)+1;
% param.config.default{cur_idx} = default_radar_params_2019_Arctic_Polar6_snow();
% param.config.base_dir{cur_idx} = '/work/ollie/ajutila/Data/mnt/HDD5/';
% param.config.config_folder_names{cur_idx} = fullfile('1903080101','UWBM');
% param.config.board_folder_names{cur_idx} = fullfile('1903080101','UWBM','%b');
% param.config.date_strs{cur_idx} = '20190308';
