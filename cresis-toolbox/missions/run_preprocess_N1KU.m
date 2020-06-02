% script run_preprocess_N1KU.m
%
% Support script for run_preprocess.m

param.config.default = [];

%% SNOW SINGLE DAY
% cur_idx = length(param.config.default)+1;
% param.config.default{cur_idx} = default_radar_params_2019_SouthDakota_N1KU_snow();
% param.config.base_dir{cur_idx} = '/cresis/snfs1/data/SnowRadar/2020_SouthDakota_CESSNA/';
% param.config.config_folder_names{cur_idx} = '20200128';
% param.config.board_folder_names{cur_idx} = '20200128';
% param.config.date_strs{cur_idx} = '20200128';

% return;

%% SNOW MULTIPLE DAYS
% date_strs = {'20200128','20200129','20200131','20200201','20200202','20200204','20200205','20200208','20200209','20200210'};
date_strs = {'20200129','20200131','20200201','20200202','20200204','20200205','20200208','20200209','20200210'};
config_format_str = '%s/';
board_format_str = '%s/';
defaults_fh = @default_radar_params_2019_SouthDakota_N1KU_snow;
base_dir = '/cresis/snfs1/data/SnowRadar/2020_SouthDakota_CESSNA/';

for idx = 1:length(date_strs)
  cur_idx = length(param.config.default)+1;
  param.config.default{cur_idx} = defaults_fh();
  param.config.base_dir{cur_idx} = base_dir;
  param.config.config_folder_names{cur_idx} = sprintf(config_format_str,date_strs{idx});
  param.config.board_folder_names{cur_idx} = sprintf(board_format_str,date_strs{idx});
  param.config.date_strs{cur_idx} = date_strs{idx};
end
