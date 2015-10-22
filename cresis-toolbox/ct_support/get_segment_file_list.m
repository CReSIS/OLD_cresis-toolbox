function [base_dir,adc_folder_name,fns,file_idxs] = get_segment_file_list(param,adc,silent_mode)
% [base_dir,adc_folder_name,fns,file_idxs] = get_segment_file_list(param,adc,silent_mode)
%
% Support function for create_vectors_* and create_records_*.
% Can also be used to get all the file information for every segment
% using run_get_segment_file_list.m.
%
% param = struct from param spreadsheet read in (read_param_xls)
% adc = optional parameter used with some radars that have multiple adcs
%   This is 1-indexed adc.
%
% Author: John Paden
%
% See also run_get_segment_file_list.m

if ~exist('silent_mode','var')
  silent_mode = false;
end

if ~isfield(param.vectors.file,'file_midfix')
  param.vectors.file.file_midfix = '';
end

if any(strcmpi(param.radar_name,{'accum','snow','kuband'}))
  adc_folder_name = param.vectors.file.adc_folder_name;
  ext = '.dat';
elseif any(strcmpi(param.radar_name,{'acords'}))
  adc_folder_name = param.vectors.file.adc_folder_name;
  ext = '';
  file_regexp = '\.[0-9]*$';
elseif strcmpi(param.radar_name,'mcrds')
  adc_folder_name = param.vectors.file.adc_folder_name;
  ext = '.raw';
elseif any(strcmpi(param.radar_name,{'accum2','snow2','kuband2','snow3','kuband3','mcords','mcords4','mcords5','snow5'}))
  % Create sub-folder name for the particular receiver
  adc_folder_name = param.vectors.file.adc_folder_name;
  board_idx_insert_idxs = strfind(adc_folder_name,'%d');
  mat_cmd = 'adc_folder_name = sprintf(adc_folder_name';
  for board_idx_insert_idx = board_idx_insert_idxs
    mat_cmd = [mat_cmd sprintf(', %d',adc)];
    adc_folder_name(board_idx_insert_idxs+1) = 'd';
  end
  mat_cmd = [mat_cmd ');'];
  eval(mat_cmd);
  
  if param.records.file_version == 401
    ext = '.dat';
  else
    ext = '.bin';
  end
elseif any(strcmpi(param.radar_name,{'mcords2','mcords3'}))
  % Create sub-folder name for the particular receiver
  adc_folder_name = param.vectors.file.adc_folder_name;
  board_idx_insert_idxs = strfind(adc_folder_name,'%b');
  mat_cmd = 'adc_folder_name = sprintf(adc_folder_name';
  board = adc_to_board(param.radar_name,adc);
  for board_idx_insert_idx = board_idx_insert_idxs
    mat_cmd = [mat_cmd sprintf(', %d',board)];
    adc_folder_name(board_idx_insert_idxs+1) = 'd';
  end
  mat_cmd = [mat_cmd ');'];
  eval(mat_cmd);
  
  ext = '.bin';
else
  error('Unsupported radar %s', param.radar_name);
end

if ~exist('file_regexp','var')
  file_regexp = '';
end

base_dir = fullfile(ct_filename_data(param,param.vectors.file.base_dir),adc_folder_name);

if nargout > 2
  if ~silent_mode
    fprintf('Getting files for %s (%s)\n', base_dir, datestr(now));
  end
  get_fns_param.regexp = file_regexp;
  fns = get_filenames(base_dir,param.vectors.file.file_prefix,param.vectors.file.file_midfix,ext,get_fns_param);
  
  % Sort ACORDS filenames because the extenions are not a standard length
  if any(strcmpi(param.radar_name,{'acords'}))
    basenames = {};
    file_idxs = [];
    new_fns = {};
    for fidx = 1:length(fns)
      fname = fname_info_acords(fns{fidx});
      new_fns{fidx} = [fname.basename sprintf('.%03d',fname.file_idx)];
    end
    [new_fns,sorted_idxs] = sort(new_fns);
    fns = fns(sorted_idxs);
  end
  
  if ~silent_mode
    fprintf('  Found %d files in %s\n', length(fns), base_dir);
  end
  
  if isempty(fns)
    fprintf('No files match the mask:\n');
    fprintf('  path: %s\n', base_dir);
    fprintf('  mask: %s*%s\n', param.vectors.file.file_prefix, ext);
    error('No files found');
  end
  
  if param.vectors.file.stop_idx == inf
    % A stop index of infinity says to include all files
    stop_idx = length(fns);
  elseif param.vectors.file.stop_idx > length(fns)
    warning('Stop index is larger than number of files available. This can be caused by an error in the stop index or missing files. dbcont to continue.');
    keyboard
    stop_idx = length(fns);
  elseif param.vectors.file.stop_idx < 0;
    % A stop index of -N says to include all but the last N files
    stop_idx = length(fns) + param.vectors.file.stop_idx;
  else
    stop_idx = param.vectors.file.stop_idx;
  end
  file_idxs = param.vectors.file.start_idx:stop_idx;
  
  if isempty(file_idxs)
    error('No files selected to load out of %i files', length(fns));
  end
end

return;
