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

%% Create each task
% =========================================================================

fun = 'preprocess_task';

ctrl = cluster_new_batch(param_override);

cluster_compile(fun,[],0,ctrl);

sparam = [];
sparam.task_function = fun;
sparam.argsin{1} = param; % Static parameters
sparam.num_args_out = 1;
for config_idx = 1:numel(param.config.default)

  cparam = merge_structs(param,param.config.default{config_idx});
  cparam.config.base_dir = param.config.base_dir{config_idx};
  cparam.config.config_folder_names = param.config.config_folder_names{config_idx};
  cparam.config.board_folder_names = param.config.board_folder_names{config_idx};
  cparam.config.date_str = param.config.date_strs{config_idx};
  cparam.config = rmfield(cparam.config,'date_strs');
  cparam.config = rmfield(cparam.config,'default');
  
  %% Input checks
  % =========================================================================
  
  if strcmpi(cparam.config.daq_type,'cresis')
    % CReSIS DAQ only parameters
    if ~isfield(cparam.config.cresis,'gps_file_mask') || isempty(cparam.config.cresis.gps_file_mask)
      % File mask relative to param.config.config_folder_name for GPS files
      % Leave empty if there are no GPS files to copy.
      cparam.config.cresis.gps_file_mask = '';
    end
  end
  
  if ~isfield(cparam.config,'date_strs') || isempty(cparam.config.date_strs)
    % Segment ID is created from the date of the config folder unless this
    % cell array is set.
    cparam.config.date_strs = {};
  end
  if ~isfield(cparam.config,'file') || isempty(cparam.config.file)
    cparam.config.file = [];
  end
  if ~isfield(cparam.config.file,'prefix') || isempty(cparam.config.file.prefix)
    cparam.config.file.prefix = '';
  end
  if ~isfield(cparam.config.file,'midfix') || isempty(cparam.config.file.midfix)
    cparam.config.file.midfix = '';
  end
  if ~isfield(cparam.config.file,'suffix') || isempty(cparam.config.file.suffix)
    cparam.config.file.suffix = '';
  end
  if ~isfield(cparam.config.file,'regexp') || isempty(cparam.config.file.regexp)
    cparam.config.file.regexp = '';
  end
  if ~isfield(cparam.config,'mat_or_bin_hdr_output') || isempty(cparam.config.mat_or_bin_hdr_output)
    cparam.config.mat_or_bin_hdr_output = '.mat';
  end
  if ~isfield(cparam.config,'max_time_gap') || isempty(cparam.config.max_time_gap)
    % Maximum time in seconds between two data records before forcing a segment break
    cparam.config.max_time_gap = 10;
  end
  if ~isfield(cparam.config,'min_seg_size') || isempty(cparam.config.min_seg_size)
    % Minimum number of files in a segment (segments with less files will be
    % discarded)
    cparam.config.min_seg_size = 2;
  end
  if ~isfield(cparam.config,'online_mode') || isempty(cparam.config.online_mode)
    cparam.config.online_mode = 0;
    % 0: not online mode
    % 1: online mode (process all files)
    % 2: online mode (process only the most recent file)
  end
  if ~isfield(cparam.config,'param_fn') || isempty(cparam.config.param_fn)
    cparam.config.param_fn ...
      = ct_filename_param(sprintf('%s_param_%s.xls',ct_output_dir(cparam.radar_name),cparam.season_name));
  end
  if ~isfield(cparam.config,'reuse_tmp_files') || isempty(cparam.config.reuse_tmp_files)
    cparam.config.reuse_tmp_files = true;
  end
  if ~isfield(cparam.config,'skip_files_end') || isempty(cparam.config.skip_files_end)
    % Number of files to ignore at the end of the segment
    cparam.config.skip_files_end = 0;
  end
  if ~isfield(cparam.config,'skip_files_start') || isempty(cparam.config.skip_files_start)
    % Number of files to ignore at the beginning of the segment
    cparam.config.skip_files_start = 0;
  end

  %% Setup preprocess task
  
  dparam = [];
  dparam.argsin{1}.season_name = cparam.season_name;
  dparam.argsin{1}.radar_name = cparam.radar_name;
  dparam.argsin{1}.config = cparam.config;
  dparam.argsin{1}.config.config_folder_name = cparam.config.config_folder_names;
  dparam.argsin{1}.config.board_folder_name = cparam.config.board_folder_names;
  dparam.argsin{1}.config.date_str = cparam.config.date_str;
  
  if strcmpi(cparam.config.daq_type,'arena')
    fns = get_filenames(fullfile(cparam.config.base_dir,cparam.config.config_folder_names),'','','.dat', ...
      struct('recursive',true,'regexp','[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]'));
    num_fns = length(fns);
    
  elseif strcmpi(cparam.config.daq_type,'cresis')
    num_fns = 0;
    for board_idx = 1:length(cparam.config.board_map)
      board = cparam.config.board_map{board_idx};
      board_folder_name = cparam.config.board_folder_names;
      board_folder_name = regexprep(board_folder_name,'%b',board);
      
      get_filenames_param = struct('regexp',cparam.config.file.regexp);
      fns = get_filenames(fullfile(cparam.config.base_dir,board_folder_name), ...
        cparam.config.file.prefix, cparam.config.file.midfix, ...
        cparam.config.file.suffix, get_filenames_param);
      num_fns = num_fns + length(fns);
    end
    
  elseif strcmpi(cparam.config.daq_type,'utua')
    num_fns = 0;
    for board_idx = 1:length(cparam.config.board_map)
      board = cparam.config.board_map{board_idx};
      board_folder_name = cparam.config.board_folder_names;
      board_folder_name = regexprep(board_folder_name,'%b',board);
      
      get_filenames_param = struct('regexp',cparam.config.file.regexp);
      fns = get_filenames(fullfile(cparam.config.base_dir,board_folder_name), ...
        cparam.config.file.prefix, cparam.config.file.midfix, ...
        cparam.config.file.suffix, get_filenames_param);
      num_fns = num_fns + length(fns);
    end
    
  else
    error('cparam.config.daq_type = %s is not a supported type', cparam.config.daq_type);
    
  end
  dparam.cpu_time = 60 + 10*num_fns;
  dparam.mem = 4e9;
  dparam.notes = sprintf('%s:%s:%s %d files', ...
    mfilename, cparam.config.base_dir, cparam.config.config_folder_names, num_fns);
  ctrl = cluster_new_task(ctrl,sparam,dparam);
end

ctrl_chain = {{ctrl}};

fprintf('Done %s\n', datestr(now));

return;
