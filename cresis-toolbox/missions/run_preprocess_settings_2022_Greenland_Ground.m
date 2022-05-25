% script run_preprocess_2022_Greenland_Ground.m
%
% Support script for run_preprocess.m
%
% Preprocess setup script for 2022_Greenland_Ground.

param.config.default = [];

%% ACCUM3 SINGLE DAY
cur_idx = length(param.config.default)+1;
param.config.default{cur_idx} = default_radar_params_2022_Greenland_Ground_accum();
param.config.base_dir{cur_idx} = '/data/';
param.config.config_folder_names{cur_idx} = '20220524';
param.config.board_folder_names{cur_idx} = '20220524/%b';
param.config.date_strs{cur_idx} = '20220524';
