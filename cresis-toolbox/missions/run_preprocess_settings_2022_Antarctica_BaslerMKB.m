% script run_preprocess_settings_2022_Antarctica_GroundGHOST.m
%
% Support script for run_preprocess.m
%
% Preprocess setup script for 2022_Antarctica_GroundGHOST.

param.config.default = [];

%% RDS SINGLE DAY
cur_idx = length(param.config.default)+1;
param.config.default{cur_idx} = @default_radar_params_2022_Antarctica_BaslerMKB_accum;
param.config.base_dir{cur_idx} = '/data/';
% param.config.config_folder_names{cur_idx} = '20221225a';
% param.config.board_folder_names{cur_idx} = '20221225a/%b';
% param.config.date_str{cur_idx} = '20221224';
param.config.config_folder_names{cur_idx} = '20230110';
param.config.board_folder_names{cur_idx} = '20230110/%b';
param.config.date_str{cur_idx} = '20230110';

% 20230109_191654_accum_digrx0_0113.dat
