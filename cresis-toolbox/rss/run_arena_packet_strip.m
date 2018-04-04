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
% See also: basic_load_arena.m

%% User Setup
% =========================================================================
param_override = [];

% Set param.radar_name and param.season_name
param = default_radar_params_2016_Greenland_TOdtu;

if ispc
  param.arena_packet_strip.base_dir = 'HF_Sounder/2016_Greenland_TO/';
else
  param.arena_packet_strip.base_dir = '/cresis/snfs1/data/HF_Sounder/2016_Greenland_TO/';
end
param.arena_packet_strip.adc_folder_names = {'20161101','20161102A','20161107A','20161107B','20161107C','20161108A','20161108B','20161110A','20161110B','20161111A','20161112A'};
param.arena_packet_strip.reuse_tmp_files = true;
param.arena_packet_strip.mat_or_bin_hdr_output = '.mat';

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









%% Automated section
% =========================================================================

fun = 'arena_packet_strip_task';

ctrl = cluster_new_batch(param_override);

cluster_compile(fun,[],0,ctrl);

sparam = [];
sparam.task_function = fun;
sparam.argsin{1} = param; % Static parameters
sparam.num_args_out = 1;
for adc_folder_idx = 1:numel(adc_folder_names)
  adc_folder_name = adc_folder_names{adc_folder_idx};
  dparam = [];
  dparam.argsin{1}.arena_packet_strip.adc_folder_name = adc_folder_name;
  fns = get_filenames(fullfile(param.arena_packet_strip.base_dir,adc_folder_name),'[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]','','.dat',struct('recursive',true));
  dparam.cpu_time = 60 + 40*length(fns);
  dparam.notes = sprintf('%s:%s:%s %d files', ...
      mfilename, param.arena_packet_strip.base_dir, adc_folder_name, length(fns));
  ctrl = cluster_new_task(ctrl,sparam,dparam);
end

ctrl_chain = {{ctrl}};

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);

% Potentially stop and inspect cluster_print_chain output to adjust
% cluster control parameters before running or to run the next lines on a
% different computer (the save/load functions are for this purpose).

return
ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.type','debug');
ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.desired_time_per_job',0);
%ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.cpu_time_mult',2);
%ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.mem_mult',2);

[ctrl_chain,chain_fn] = cluster_load_chain([],chain_id);
ctrl_chain = cluster_run(ctrl_chain);
