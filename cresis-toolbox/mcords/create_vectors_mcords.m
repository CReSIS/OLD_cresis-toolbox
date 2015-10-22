function create_vectors_mcords(param,param_override)
% create_vectors_mcords(param,param_override)
%
% Creates a .mat file containg a structure which has a time stamp,
% position, and filename for every filename passed in. It is used
% with the plot_vectors command.
%
% Can be run as a function by passing in the param argument
% or a script (by setting the default value of param).
%
% Output file contains:
% vectors: structure with the following fields
%    .year: vector
%    .month: vector
%    .day: vector
%    .timeUTC: vector of seconds since epoch Jan 1, 1970
%    .lat: latitude (double vector, degrees)
%    .lon: longitude (double vector, degrees)
%    .elevation: elevation (double vector, meters)
%    .file: structure containing fields to pass to load_MCoRDS.m
%
% Inputs:
% param = struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Examples:
%   See default user settings inside this function.
%
% Authors: John Paden
%
% See also: create_vectors_snow.m, plot_vectors.m, master.m, make_gps_*.m

% =====================================================================
% General Setup
% =====================================================================

dbstack_info = dbstack;
if ~exist('param','var') || isempty(param) || length(dbstack_info) == 1
  % =====================================================================
  % Debug Setup
  % =====================================================================
  param = read_param_xls(ct_filename_param('rds_param_2011_Greenland_P3.xls'),'20110310_03');
  
  clear('param_override');
  param_override.sched.type = 'no scheduler';
  param_override.sched.rerun_only = true;

  % Input checking
  if ~exist('param','var')
    error('A struct array of parameters must be passed in\n');
  end
  global gRadar;
  if exist('param_override','var')
    param_override = merge_structs(gRadar,param_override);
  else
    param_override = gRadar;
  end
  
elseif ~isstruct(param)
  % Functional form
  param();
end
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

% =====================================================================
% Get the files
% =====================================================================

[base_dir,adc_folder_name,fns,file_idxs] = get_segment_file_list(param,param.vectors.file.adc);

clear vectors;
vectors_idx = 0;

% Read header information
file_param.clk = param.radar.fs;
for file_num = file_idxs
  fn = fns{file_num};
  
  % Find first good header
  fname = fname_info_mcords(fn);
  file_param.clk = param.radar.fs;
  
  try
    if fname.file_idx == 0
      % First file after record starts has junk data at start:
      % Skip first 64 MB (should just be 32 MB, but a lot of errors
      % seem to occur in the first 64 MB)
      finfo = frame_sync_info(fn,struct('fast',1,'first_byte',2^26));
    else
      finfo = frame_sync_info(fn,struct('fast',1));
    end
  catch
    % Sync failed (probably too short of a file)
    continue
  end;
  data_start = finfo.ideal_sync_pos(1);
  
  % Read first header
  [fid msg] = fopen(fn,'r','ieee-be');
  if fid == -1
    fprintf('Unable to open file: %s\n',fn);
    error('File:Open',msg);
  end
  fseek(fid,data_start+8,'bof');
  
  vectors_idx = vectors_idx + 1;
  hdr_time_sod(vectors_idx) = double(fread(fid,1,'uint32')) + double(fread(fid,1,'uint32'))/(param.radar.fs/2);
  
  fclose(fid);
  
  % Store the filename information
  fprintf('  File index %d (filename idx %d), (%s)\n', ...
    file_num, fname.file_idx, datestr(now,'HH:MM:SS'));
  vectors.file(vectors_idx).idx = fname.file_idx;
  vectors.file(vectors_idx).rxs = fname.rx;
  vectors.file(vectors_idx).recording_group = fname.group;
  vectors.file(vectors_idx).base_dir = base_dir;
  vectors.file(vectors_idx).adc_folder_name = adc_folder_name;
  vectors.file(vectors_idx).datenum = fname.datenum;
  vectors.fileNumber(vectors_idx) = file_num;
  
end

vectors = sync_radar_to_gps(param,vectors,hdr_time_sod);

% Make sure output directory exists
[out_dir out_name] = fileparts(ct_filename_support(param,param.vectors.out_fn,'vectors'));
if ~exist(out_dir,'dir')
  mkdir(out_dir);
end

param_vectors = param;
vectors_fn = ct_filename_support(param,param.vectors.out_fn,'vectors');
fprintf('  Saving file %s\n', vectors_fn);
save(vectors_fn,'vectors','param_vectors');

return;
