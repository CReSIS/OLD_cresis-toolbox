function ctrl_chain = preprocess(param,param_override)
% ctrl_chain = preprocess(param,param_override)
%
% 1. Copy metadata files if applicable
% 2. Load raw data files:
%   a. Create temporary header files for create_records
%   b. Create network header stripped raw data files if applicable
% 3. Break data into segments
% 4. Create parameter spreadsheet entries
% 5. Load GPS data and create maps of each segment and file if applicable
%
% Author: John Paden
%
% See also: run_preprocess.m, preprocess.m, preprocess_task.m

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, '', datestr(now));
fprintf('=====================================================================\n');

%% Input checks
% =========================================================================

if ~isfield(param.config,'file') || isempty(param.config.file)
  param.config.file = [];
end

if ~isfield(param.config.file,'prefix') || isempty(param.config.file.prefix)
  param.config.file.prefix = '';
end
if ~isfield(param.config.file,'midfix') || isempty(param.config.file.midfix)
  param.config.file.midfix = '';
end
if ~isfield(param.config.file,'suffix') || isempty(param.config.file.suffix)
  param.config.file.suffix = '';
end
if ~isfield(param.config.file,'regexp') || isempty(param.config.file.regexp)
  param.config.file.regexp = '';
end
if ~isfield(param.config,'online_mode') || isempty(param.config.online_mode)
  param.config.online_mode = 0;
  % 0: not online mode
  % 1: online mode (process all files)
  % 2: online mode (process only the most recent file)
end

if ~isfield(param.config,'max_time_gap') || isempty(param.config.max_time_gap)
  % Maximum time in seconds between two data records before forcing a segment break
  param.config.max_time_gap = 10;
end
if ~isfield(param.config,'min_seg_size') || isempty(param.config.min_seg_size)
  % Minimum number of files in a segment (segments with less files will be
  % discarded)
  param.config.min_seg_size = 2;
end
if ~isfield(param.config,'skip_files_start') || isempty(param.config.skip_files_start)
  % Number of files to ignore at the beginning of the segment
  param.config.skip_files_start = 0;
end
if ~isfield(param.config,'skip_files_end') || isempty(param.config.skip_files_end)
  % Number of files to ignore at the end of the segment
  param.config.skip_files_end = 0;
end
if ~isfield(param.config,'date_strs') || isempty(param.config.date_strs)
  % Segment ID is created from the date of the config folder unless this
  % cell array is set.
  param.config.date_strs = {};
end
if ~isfield(param.config,'reuse_tmp_files') || isempty(param.config.reuse_tmp_files)
  param.config.reuse_tmp_files = true;
end
if ~isfield(param.config,'mat_or_bin_hdr_output') || isempty(param.config.mat_or_bin_hdr_output)
  param.config.mat_or_bin_hdr_output = '.mat';
end
if ~isfield(param.config,'param_fn') || isempty(param.config.param_fn)
  param.config.param_fn ...
    = ct_filename_param(sprintf('%s_param_%s.xls',ct_output_dir(param.radar_name),param.season_name));
end


%% Create each task
% =========================================================================

fun = 'preprocess_task';

ctrl = cluster_new_batch(param_override);

cluster_compile(fun,[],0,ctrl);

sparam = [];
sparam.task_function = fun;
sparam.argsin{1} = param; % Static parameters
sparam.num_args_out = 1;
for config_folder_idx = 1:numel(param.config.config_folder_names)
  config_folder_name = param.config.config_folder_names{config_folder_idx};
  board_folder_name = param.config.board_folder_names{config_folder_idx};
  
  dparam = [];
  dparam.argsin{1}.config = param.config;
  dparam.argsin{1}.config.config_folder_name = config_folder_name;
  dparam.argsin{1}.config.board_folder_name = board_folder_name;
  if numel(param.config.date_strs) >= config_folder_idx
    dparam.argsin{1}.config.date_str = param.config.date_strs{config_folder_idx};
  end
  
  if strcmpi(param.config.daq_type,'arena')
    fns = get_filenames(fullfile(param.config.base_dir,config_folder_name),'','','.dat', ...
      struct('recursive',true,'regexp','[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]'));
    num_fns = length(fns);
    
  elseif strcmpi(param.config.daq_type,'cresis')
    num_fns = 0;
    for board_idx = 1:length(param.config.daq.board_map)
      board = param.config.daq.board_map{board_idx};
      board_folder_name = param.config.board_folder_names{config_folder_idx};
      board_folder_name = regexprep(board_folder_name,'%b',board);
      
      get_filenames_param = struct('regexp',param.config.file.regexp);
      fns = get_filenames(fullfile(param.config.base_dir,board_folder_name), ...
        param.config.file.prefix, param.config.file.midfix, ...
        param.config.file.suffix, get_filenames_param);
      num_fns = num_fns + length(fns);
    end
    
  else
    error('param.config.daq_type = %s is not a supported type', param.config.daq_type);
    
  end
  dparam.cpu_time = 60 + 10*num_fns;
  dparam.notes = sprintf('%s:%s:%s %d files', ...
    mfilename, param.config.base_dir, config_folder_name, num_fns);
  ctrl = cluster_new_task(ctrl,sparam,dparam);
end

ctrl_chain = {{ctrl}};

fprintf('Done %s\n', datestr(now));

return;
