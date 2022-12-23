% script run_preprocess_2022_Antarctica_GroundGHOST.m
%
% Support script for run_preprocess.m
%
% Preprocess setup script for 2022_Antarctica_GroundGHOST.

param.config.default = [];

%% RDS SINGLE DAY
cur_idx = length(param.config.default)+1;
param.config.default{cur_idx} = @default_radar_params_2022_Antarctica_GroundGHOST_rds;
param.config.base_dir{cur_idx} = '/data/';
param.config.config_folder_names{cur_idx} = '20221223';
param.config.board_folder_names{cur_idx} = '20221223/%b';
param.config.date_str{cur_idx} = '20221223';
