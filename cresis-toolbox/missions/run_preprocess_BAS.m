% script run_preprocess_BAS.m
%
% Support script for run_preprocess.m

param.config.default = [];

%% ACCUM3 TObas
cur_idx = length(param.config.default)+1;
param.config.default{cur_idx} = default_radar_params_2018_Antarctica_TObas_accum();
param.config.base_dir{cur_idx} = '/data/';
param.config.config_folder_names{cur_idx} = '20181004';
param.config.board_folder_names{cur_idx} = '20181004/%b';
param.config.date_strs{cur_idx} = '20181004';
