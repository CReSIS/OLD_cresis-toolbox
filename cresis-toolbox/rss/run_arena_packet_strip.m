% script run_arena_packet_strip.m
%
% Strips data out of Arena packets and creates:
% 1. Files of just radar data
% 2. Files with just header data
%   a. binary flat file .hdr
%   b. Matlab .mat format (if supported)
% NOTE: This script only works on files with fixed header lengths. More
% specifically, every header must have the payload length field in the same spot.
%
% Author: John Paden
%
% See also: basic_load_arena.m, run_arena_packet_strip.m,
% arena_packet_strip.m, arena_packet_strip_task.m

%% User Setup
% =========================================================================
param_override = [];

% Set param.radar_name and param.season_name and get radar default
% parameters.
% [param,defaults] = default_radar_params_2018_Antarctica_TObas_accum3;
[param,defaults] = default_radar_params_2018_Antarctica_Ground_mcords6;

param.arena_packet_strip.defaults = defaults;
param.arena_packet_strip.default_param = param;

if ispc
  param.arena_packet_strip.base_dir = 'E:\tmp\2018_Antarctica_TObas\';
else
  param.arena_packet_strip.base_dir = '/mnt/scratch/';
end
% param.arena_packet_strip.config_folder_names = {'20180817'};
% param.arena_packet_strip.board_folder_names = {'20180817/%b'};
param.arena_packet_strip.config_folder_names = {'20180829'};
param.arena_packet_strip.board_folder_names = {'20180829/%b'};
param.arena_packet_strip.board_map = {'digrx0','digrx1','digrx2','digrx3'};
param.arena_packet_strip.tx_map = {'awg0','awg1','awg2','awg3'};
param.arena_packet_strip.reuse_tmp_files = true;
param.arena_packet_strip.mat_or_bin_hdr_output = '.mat';
% param.arena_packet_strip.param_fn = ct_filename_param('accum_param_2018_Antarctica_TObas.xls');
param.arena_packet_strip.param_fn = ct_filename_param('rds_param_2018_Antarctica_Ground.xls');

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

ctrl_chain = arena_packet_strip(param,param_override);

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
