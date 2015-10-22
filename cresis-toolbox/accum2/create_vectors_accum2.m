function create_vectors_accum2(param,param_override)
% create_vectors_accum2(param,param_override)
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
%    .file: structure containing fields to pass to load_accum.m
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
% Authors: John Paden, Logan Smith
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
  param = read_param_xls(ct_filename_param('accum_param_2011_Greenland_P3.xls'),'20110310_03');
  
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

[base_dir,adc_folder_name,fns,file_idxs] = get_segment_file_list(param);

clear vectors;
vectors_idx = 0;
  file_param.clk = param.radar.fs;
for file_num = 1:length(file_idxs)
  fn = fns{file_idxs(file_num)};
  fname = fname_info_accum2(fn);
  
  fprintf('  File index %d/%d (filename idx %d) (%s)\n', ...
    file_num, length(file_idxs), fname.file_idx, datestr(now,'HH:MM:SS'));

  
  % Read header information
  hdr = basic_load_accum2(fn, struct('clk',param.radar.fs));
  
  vectors_idx = vectors_idx + 1;
  % Store the filename information
  vectors.file(vectors_idx).idx = fname.file_idx;
  vectors.file(vectors_idx).adcs = 1;
  vectors.file(vectors_idx).recording_group = fname.group;
  vectors.file(vectors_idx).base_dir = param.vectors.file.base_dir;
  vectors.file(vectors_idx).adc_folder_name = param.vectors.file.adc_folder_name;
  vectors.file(vectors_idx).datenum = fname.datenum;
  vectors.fileNumber(vectors_idx) = file_num;
  vectors.comp_time(vectors_idx) = hdr.comp_time;
  vectors.radar_time(vectors_idx) = hdr.radar_time;
end

%% Synchronize to GPS data
vectors = sync_radar_to_gps(param,vectors,vectors.radar_time,vectors.comp_time);

% Make sure output directory exists
% Make sure output directory exists
out_fn = ct_filename_support(param,param.vectors.out_fn,'vectors');
[out_dir out_name] = fileparts(out_fn);
if ~exist(out_dir,'dir')
  mkdir(out_dir);
end

param_vectors = param.vectors;

fprintf('  Saving file %s\n', out_fn);
save(out_fn,'vectors','param_vectors');

return;










