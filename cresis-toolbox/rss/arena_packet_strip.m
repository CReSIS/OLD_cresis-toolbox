function ctrl_chain = arena_packet_strip(param,param_override)
% ctrl_chain = arena_packet_strip(param,param_override)
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

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, '', datestr(now));
fprintf('=====================================================================\n');

%% Automated section
% =========================================================================

fun = 'arena_packet_strip_task';

ctrl = cluster_new_batch(param_override);

cluster_compile(fun,[],0,ctrl);

sparam = [];
sparam.task_function = fun;
sparam.argsin{1} = param; % Static parameters
sparam.num_args_out = 1;
for config_folder_idx = 1:numel(param.arena_packet_strip.config_folder_names)
  config_folder_name = param.arena_packet_strip.config_folder_names{config_folder_idx};
  board_folder_name = param.arena_packet_strip.board_folder_names{config_folder_idx};
  dparam = [];
  dparam.argsin{1}.arena_packet_strip.config_folder_name = config_folder_name;
  dparam.argsin{1}.arena_packet_strip.board_folder_name = board_folder_name;
  
  if ~isfield(param.arena_packet_strip,'board_map') || isempty(param.arena_packet_strip.board_map)
    dparam.argsin{1}.arena_packet_strip.board_map = param.arena_packet_strip.defaults{1}.board_map;
  else
    dparam.argsin{1}.arena_packet_strip.board_map = param.arena_packet_strip.board_map;
  end
  if ~isfield(param.arena_packet_strip,'tx_map') || isempty(param.arena_packet_strip.tx_map)
    dparam.argsin{1}.arena_packet_strip.tx_map = param.arena_packet_strip.defaults{1}.tx_map;
  else
    dparam.argsin{1}.arena_packet_strip.tx_map = param.arena_packet_strip.tx_map;
  end
  if ~isfield(param.arena_packet_strip,'reuse_tmp_files') || isempty(param.arena_packet_strip.reuse_tmp_files)
    dparam.argsin{1}.arena_packet_strip.reuse_tmp_files = true;
  end
  if ~isfield(param.arena_packet_strip,'mat_or_bin_hdr_output') || isempty(param.arena_packet_strip.mat_or_bin_hdr_output)
    dparam.argsin{1}.arena_packet_strip.mat_or_bin_hdr_output = '.mat';
  end
  if ~isfield(param.arena_packet_strip,'param_fn') || isempty(param.arena_packet_strip.param_fn)
    dparam.argsin{1}.arena_packet_strip.param_fn ...
      = ct_filename_param(sprintf('%s_param_%s.xls',ct_output_dir(param.radar_name),param.season_name));
  end
  
  fns = get_filenames(fullfile(param.arena_packet_strip.base_dir,config_folder_name),'','','.dat',struct('recursive',true,'regexp','[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]'));
  dparam.cpu_time = 60 + 40*length(fns);
  dparam.notes = sprintf('%s:%s:%s %d files', ...
      mfilename, param.arena_packet_strip.base_dir, config_folder_name, length(fns));
  ctrl = cluster_new_task(ctrl,sparam,dparam);
end

ctrl_chain = {{ctrl}};
    
fprintf('Done %s\n', datestr(now));

return;
