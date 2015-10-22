function create_vectors_fmcw(param,param_override)
% create_vectors_fmcw(param,param_override)
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
  param = read_param_xls(ct_filename_param('kuband_param_2013_Greenland_P3.xls'),'20130315_01');
  
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
for file_num = 1:length(file_idxs)
  fn = fns{file_idxs(file_num)};
  fname = fname_info_fmcw(fn);
  
  fprintf('  File index %d/%d (filename idx %d) (%s)\n', ...
    file_num, length(file_idxs), fname.file_idx, datestr(now,'HH:MM:SS'));
  
  % Read header information
  if any(strcmpi(param.radar_name,{'snow','kuband'}))
    hdr = basic_load_fmcw(fn, struct('clk',param.radar.fs, ...
      'utc_time_halved',param.vectors.gps.utc_time_halved));
  elseif any(strcmpi(param.radar_name,{'snow2','kuband2'}))
    hdr = basic_load_fmcw2(fn, struct('clk',param.radar.fs,'file_version', param.records.file_version));
  elseif any(strcmpi(param.radar_name,{'snow3','kuband3'}))
    if param.records.file_version == 4
      hdr = basic_load_fmcw2(fn, struct('clk',param.radar.fs,'file_version', param.records.file_version));
    else
      hdr = basic_load_fmcw3(fn, struct('clk',param.radar.fs));
    end
  end
  
  vectors_idx = vectors_idx + 1;
  % Store the filename information
  vectors.file(vectors_idx).idx = fname.file_idx;
  vectors.file(vectors_idx).adcs = 1;
  vectors.file(vectors_idx).recording_group = fname.group;
  vectors.file(vectors_idx).base_dir = param.vectors.file.base_dir;
  vectors.file(vectors_idx).adc_folder_name = param.vectors.file.adc_folder_name;
  vectors.file(vectors_idx).datenum = fname.datenum;
  vectors.fileNumber(vectors_idx) = file_num;
  if strcmpi(param.season_name,'2014_Greenland_P3') & (strcmpi(param.day_seg,'20140421_02') |...
      strcmpi(param.day_seg,'20140423_01') | strcmpi(param.day_seg,'20140502_00') |...
      strcmpi(param.day_seg,'20140508_02')) 
    hdr_time_sod(vectors_idx) = str2num(fn(end-14:end-13))*3600 + str2num(fn(end-12:end-11))*60 + str2num(fn(end-10:end-9));
  else
    hdr_time_sod(vectors_idx) = hdr.utc_time_sod;
  end
end
if strcmpi(param.season_name,'2015_Greenland_C130')
  if strcmpi(param.day_seg,'20150324_04') 
    hdr_time_sod(7) = hdr_time_sod(6) + (hdr_time_sod(6) - hdr_time_sod(5));
  end
  if strcmpi(param.day_seg,'20150410_01')
    hdr_time_sod(405) = hdr_time_sod(404) + (hdr_time_sod(404) - hdr_time_sod(403));
    hdr_time_sod(591) = hdr_time_sod(590) + (hdr_time_sod(590) - hdr_time_sod(589));
    hdr_time_sod(593) = hdr_time_sod(592) + (hdr_time_sod(590) - hdr_time_sod(589));
  end
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
