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

if ~isfield(param.preprocess,'file') || isempty(param.preprocess.file)
  param.preprocess.file = [];
end

if ~isfield(param.preprocess.file,'prefix') || isempty(param.preprocess.file.prefix)
  param.preprocess.file.prefix = '';
end
if ~isfield(param.preprocess.file,'midfix') || isempty(param.preprocess.file.midfix)
  param.preprocess.file.midfix = '';
end
if ~isfield(param.preprocess.file,'suffix') || isempty(param.preprocess.file.suffix)
  param.preprocess.file.suffix = '';
end
if ~isfield(param.preprocess.file,'regexp') || isempty(param.preprocess.file.regexp)
  param.preprocess.file.regexp = '';
end
if ~isfield(param.preprocess,'online_mode') || isempty(param.preprocess.online_mode)
  param.preprocess.online_mode = 0;
  % 0: not online mode
  % 1: online mode (process all files)
  % 2: online mode (process only the most recent file)
end

if ~isfield(param.preprocess,'max_time_gap') || isempty(param.preprocess.max_time_gap)
  % Maximum time in seconds between two data records before forcing a segment break
  param.preprocess.max_time_gap = 10;
end
if ~isfield(param.preprocess,'min_seg_size') || isempty(param.preprocess.min_seg_size)
  % Minimum number of files in a segment (segments with less files will be
  % discarded)
  param.preprocess.min_seg_size = 2;
end
if ~isfield(param.preprocess,'skip_files_start') || isempty(param.preprocess.skip_files_start)
  % Number of files to ignore at the beginning of the segment
  param.preprocess.skip_files_start = 0;
end
if ~isfield(param.preprocess,'skip_files_end') || isempty(param.preprocess.skip_files_end)
  % Number of files to ignore at the end of the segment
  param.preprocess.skip_files_end = 0;
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
for config_folder_idx = 1:numel(param.preprocess.config_folder_names)
  config_folder_name = param.preprocess.config_folder_names{config_folder_idx};
  board_folder_name = param.preprocess.board_folder_names{config_folder_idx};
  
  dparam = [];
  dparam.argsin{1}.preprocess.config_folder_name = config_folder_name;
  dparam.argsin{1}.preprocess.board_folder_name = board_folder_name;
  dparam.argsin{1}.preprocess.date_str = param.preprocess.date_strs{config_folder_idx};
  dparam.argsin{1}.preprocess.daq = param.preprocess.daq;
  dparam.argsin{1}.preprocess.wg = param.preprocess.wg;
  
  if ~isfield(param.preprocess,'reuse_tmp_files') || isempty(param.preprocess.reuse_tmp_files)
    dparam.argsin{1}.preprocess.reuse_tmp_files = true;
  end
  if ~isfield(param.preprocess,'mat_or_bin_hdr_output') || isempty(param.preprocess.mat_or_bin_hdr_output)
    dparam.argsin{1}.preprocess.mat_or_bin_hdr_output = '.mat';
  end
  if ~isfield(param.preprocess,'param_fn') || isempty(param.preprocess.param_fn)
    dparam.argsin{1}.preprocess.param_fn ...
      = ct_filename_param(sprintf('%s_param_%s.xls',ct_output_dir(param.radar_name),param.season_name));
  end
  
  if strcmpi(param.preprocess.daq.type,'arena')
    fns = get_filenames(fullfile(param.preprocess.base_dir,config_folder_name),'','','.dat', ...
      struct('recursive',true,'regexp','[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]'));
    num_fns = length(fns);
    
  elseif strcmpi(param.preprocess.daq.type,'cresis')
    num_fns = 0;
    for board_idx = 1:length(param.preprocess.daq.board_map)
      board = param.preprocess.daq.board_map{board_idx};
      board_folder_name = param.preprocess.board_folder_names{config_folder_idx};
      board_folder_name = regexprep(board_folder_name,'%b',board);
      
      get_filenames_param = struct('regexp',param.preprocess.file.regexp);
      fns = get_filenames(fullfile(param.preprocess.base_dir,board_folder_name), ...
        param.preprocess.file.prefix, param.preprocess.file.midfix, ...
        param.preprocess.file.suffix, get_filenames_param);
      num_fns = num_fns + length(fns);
    end
    
  else
    error('param.preprocess.daq.type = %s is not a supported type', param.preprocess.daq.type);
    
  end
  dparam.cpu_time = 60 + 10*num_fns;
  dparam.notes = sprintf('%s:%s:%s %d files', ...
    mfilename, param.preprocess.base_dir, config_folder_name, num_fns);
  ctrl = cluster_new_task(ctrl,sparam,dparam);
end

ctrl_chain = {{ctrl}};

fprintf('Done %s\n', datestr(now));

return;
