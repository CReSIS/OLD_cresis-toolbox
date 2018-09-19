function [base_dir,board_folder_name,fns,file_idxs] = get_segment_file_list(param,board_idx)
% [base_dir,board_folder_name,fns,file_idxs] = get_segment_file_list(param,board_idx)
%
% Support function for create_vectors.m and create_records.m.
% Can also be used to get all the file information for every segment
% using run_get_segment_file_list.m.
%
% param: struct from param spreadsheet read in by read_param_xls
%  .records
%   .file
%    .version: raw file version (see "raw file guide" on wiki)
%    .board_folder_name: board folder name (%b in the filename will be
%      replaced by the board number)
%    .base_dir: base directory
%    .prefix: beginning filename search term
%    .midfix: middle filename search term
%    .regexp: additional regular expression to run after the initial file
%      search
% board: optional parameter used with some radars that have multiple adcs
%
% Author: John Paden
%
% See also: get_frame_id, get_raw_files.m, get_segment_file_list.m,
%   run_get_segment_file_list.m

if ~isfield(param.records.file,'file_midfix')
  param.records.file.midfix = '';
end

% Determine default regular expression to apply to file search
if ~isfield(param.records.file,'regexp') || isempty(param.records.file.regexp)
  if any(param.records.file.version == [405 406])
    param.records.file.regexp = '\.[0-9]*$';
  elseif any(param.records.file.version == [409])
    param.records.file.regexp = '\.[0-9][0-9][0-9]$';
  else
    param.records.file.regexp = '';
  end
end

% Determine raw filename extension and check file version
if any(param.records.file.version == [1 9:10 101 103 401 412])
  ext = '.dat';
elseif any(param.records.file.version == [2:8 102 402:408 411])
  ext = '.bin';
elseif any(param.records.file.version == [410])
  ext = '.raw';
elseif any(param.records.file.version == [405 406 409])
  ext = '';
else 
  error('Unsupported file version\n');
end

board_folder_name = param.records.file.board_folder_name;
board_folder_name = regexprep(board_folder_name,'%b',param.records.file.boards{board_idx});

base_dir = fullfile(ct_filename_data(param,param.records.file.base_dir),board_folder_name);

%% Return file list
if nargout > 2
  if 0
    fprintf('Getting files for %s (%s)\n', base_dir, datestr(now));
  end
  get_fns_param = struct('regexp',param.records.file.regexp);
  if ~isfield(param.records.file,'file_suffix')
    fns = get_filenames(base_dir,param.records.file.prefix,param.records.file.midfix,ext,get_fns_param);
  else
    fns = get_filenames(base_dir,param.records.file.prefix,param.records.file.midfix,param.records.file.suffix,get_fns_param);
  end
  
  % Sort ACORDS filenames because the extenions are not a standard length
  if any(param.records.file.version == [405 406])
    basenames = {};
    file_idxs = [];
    new_fns = {};
    for fidx = 1:length(fns)
      fname = fname_info_acords(fns{fidx},struct('hnum',1,'file_version',param.records.file.version));
      new_fns{fidx} = [fname.basename sprintf('.%03d',fname.file_idx)];
    end
    [new_fns,sorted_idxs] = sort(new_fns);
    fns = fns(sorted_idxs);
  end

  if 0
    fprintf('  Found %d files in %s\n', length(fns), base_dir);
  end
  
  if isempty(fns)
    fprintf('No files match the mask:\n');
    fprintf('  path: %s\n', base_dir);
    fprintf('  mask: %s*%s\n', param.records.file.prefix, ext);
    error('No files found');
  end
  
  if param.records.file.stop_idx == inf
    % A stop index of infinity says to include all files
    stop_idx = length(fns);
  elseif param.records.file.stop_idx > length(fns)
    warning('Stop index (%d) is larger than number of files available (%d). This can be caused by an error in the stop index or missing files. dbcont to continue.',param.records.file.stop_idx,length(fns));
    keyboard
    stop_idx = length(fns);
  elseif param.records.file.stop_idx < 0;
    % A stop index of -N says to include all but the last N files
    stop_idx = length(fns) + param.records.file.stop_idx;
  else
    stop_idx = param.records.file.stop_idx(board_idx);
  end
  file_idxs = param.records.file.start_idx(board_idx):stop_idx;
  
  if isempty(file_idxs)
    error('No files selected to load out of %i files', length(fns));
  end
end

return;
