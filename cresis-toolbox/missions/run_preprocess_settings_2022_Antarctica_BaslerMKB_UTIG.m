% script run_preprocess_settings_2022_Antarctica_GroundGHOST.m
%
% Support script for run_preprocess.m
%
% Preprocess setup script for 2022_Antarctica_GroundGHOST.

param.config.default = [];

%% RDS SINGLE DAY
cur_idx = length(param.config.default)+1;
param.config.default{cur_idx} = @default_radar_params_2022_Antarctica_BaslerMKB_rds;
param.config.base_dir{cur_idx} = '/data/UTIG/orig/xped/CXA1/acqn/MARFA/';
param.config.config_folder_names{cur_idx} = '../../ELSA/F13';
param.config.board_folder_names{cur_idx} = 'F13';
param.config.date_str{cur_idx} = '20230116';

%   fn = '/data/UTIG/orig/xped/CXA1/acqn/MARFA/F13/radar0_20230116-200145-0001.dat';
