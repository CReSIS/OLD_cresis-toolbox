% script run_preprocess.m
%
% Runs script preprocess.m
%
% Instructions:
% 1. Set your default radar parameters file
% 2. Set the input directory of the data files
% 3. Set the input directory of the config files if applicable
% 4. Run the script
%
% Author: John Paden

%% User Setup
% =========================================================================
param_override = [];

% Set param.radar_name and param.season_name and get radar default
% parameters.
% param = default_radar_params_2018_Greenland_P3_rds;
param = default_radar_params_2018_Greenland_P3_snow;
% param = default_radar_params_2018_Antarctica_TObas_accum3;
% param = default_radar_params_2018_Antarctica_Ground_mcords6;

if ispc
  param.preprocess.base_dir = 'E:\tmp\2018_Antarctica_TObas\';
else
  param.preprocess.base_dir = '/N/dcwan/projects/cresis/2018_Greenland_P3/';
  % param.preprocess.base_dir = '/data/';
end
% param.preprocess.config_folder_names = {'20180405/mcords/'};
% param.preprocess.board_folder_names = {'20180405/mcords/%b'};
% param.preprocess.date_strs = {'20180405'};
param.preprocess.config_folder_names = {'20180405/fmcw/snow/'};
param.preprocess.board_folder_names = {'20180405/fmcw/snow/'};
param.preprocess.date_strs = {'20180405'};
% param.preprocess.config_folder_names = {'20180929/'};
% param.preprocess.board_folder_names = {'20180929/%b'};


dbstop if error;
param_override.cluster.type = 'debug';
% param_override.cluster.type = 'matlab';
% param_override.cluster.type = 'debug';
%param_override.cluster.rerun_only = true;
param_override.cluster.desired_time_per_job  = 0*60;
% param_override.cluster.cpu_time_mult  = 2;
% param_override.cluster.mem_mult  = 2;

%% Automated Section
% =====================================================================

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

ctrl_chain = preprocess(param,param_override);

cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);

% Potentially stop and inspect cluster_print_chain output to adjust
% cluster control parameters before running or to run the next lines on a
% different computer (the save/load functions are for this purpose).

return
%ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.desired_time_per_job',5*60);
%ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.cpu_time_mult',2);
%ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.mem_mult',2);

[ctrl_chain,chain_fn] = cluster_load_chain([],chain_id);
ctrl_chain = cluster_run(ctrl_chain);

