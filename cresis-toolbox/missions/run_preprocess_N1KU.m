% script run_preprocess_N1KU.m
%
% Support script for run_preprocess.m

param.config.default = [];

%% SNOW SINGLE DAY
% cur_idx = length(param.config.default)+1;
% param.config.default{cur_idx} = default_radar_params_2020_SouthDakota_N1KU_snow();
% param.config.base_dir{cur_idx} = '/cresis/snfs1/data/SnowRadar/2020_SouthDakota_CESSNA/';
% param.config.config_folder_names{cur_idx} = '20200128';
% param.config.board_folder_names{cur_idx} = '20200128';
% param.config.date_strs{cur_idx} = '20200128';

% return;

%% SNOW MULTIPLE DAYS
date_strs = {'20210201'};
config_format_str = '';
board_format_str = '';
defaults_fh = @default_radar_params_2020_SouthDakota_N1KU_snow;
base_dir = '/cresis/snfs1/projects/Field_Experiments/2021_Winter_Cessna172/Lab_tests/recording_tests/loopback_wGPS/seg2/';
base_dir = '/cresis/snfs1/projects/Field_Experiments/2021_Winter_Cessna172/Lab_tests/recording_tests/loopback_wGPS/seg5/';

for idx = 1:length(date_strs)
  cur_idx = length(param.config.default)+1;
  param.config.default{cur_idx} = defaults_fh();
  param.config.base_dir{cur_idx} = base_dir;
  param.config.config_folder_names{cur_idx} = sprintf(config_format_str,date_strs{idx});
  param.config.board_folder_names{cur_idx} = sprintf(board_format_str,date_strs{idx});
  param.config.date_strs{cur_idx} = date_strs{idx};
end
