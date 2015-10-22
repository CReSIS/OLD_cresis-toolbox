function create_vectors_accum(param,param_override)
% create_vectors_accum(param,param_override)
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
%    .file: structure containing fields to pass
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
%   See master
%
% Authors: John Paden
%
% See also: master

% =====================================================================
% General Setup
% =====================================================================

dbstack_info = dbstack;
if ~exist('param','var') || isempty(param) || length(dbstack_info) == 1
  % =====================================================================
  % Debug Setup
  % =====================================================================
  param = read_param_xls(ct_filename_param('accum_param_2011_Greenland_P3.xls'),'20110411_02');
  
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

base_dir = fullfile(param.vectors.file.base_dir,param.vectors.file.adc_folder_name);

fprintf('Getting files for %s (%s)\n', base_dir, datestr(now));
fns = get_filenames(base_dir,param.vectors.file.file_prefix,'','.dat');
fprintf('  Found %d files in %s\n', length(fns), base_dir);
if isempty(fns)
  return;
end

if param.vectors.file.stop_idx > length(fns)
  stop_idx = length(fns);
else
  stop_idx = param.vectors.file.stop_idx;
end
file_idxs = param.vectors.file.start_idx:stop_idx;

vectors_idx = 0;
for file_num = 1:length(file_idxs)
  fn = fns{file_idxs(file_num)};
  fname = fname_info_accum(fn);
  
  fprintf('  File index %d (filename idx %d) (%s)\n', ...
    file_num, fname.file_idx, datestr(now,'HH:MM:SS'));
  
  % Read header information
  hdr = basic_load_accum(fn, struct('clk',param.radar.fs));

  vectors_idx = vectors_idx + 1;
  % Store the filename information
  vectors.file(vectors_idx).idx = fname.file_idx;
  vectors.file(vectors_idx).adcs = 1;
  vectors.file(vectors_idx).recording_group = 1;
  vectors.file(vectors_idx).base_dir = param.vectors.file.base_dir;
  vectors.file(vectors_idx).adc_folder_name = param.vectors.file.adc_folder_name;
  vectors.file(vectors_idx).datenum = fname.datenum;
  vectors.fileNumber(vectors_idx) = file_num;
  hdr_time_sod(vectors_idx) = hdr.utc_time_sod;
end

vectors = sync_radar_to_gps(param,vectors,hdr_time_sod);

% Make sure output directory exists
out_fn = ct_filename_support(param,param.vectors.out_fn,'vectors');
[out_dir out_name] = fileparts(out_fn);
if ~exist(out_dir,'dir')
  mkdir(out_dir);
end

param_vectors = param;

fprintf('  Saving file %s\n', out_fn);
save(out_fn,'vectors','param_vectors');

return;
