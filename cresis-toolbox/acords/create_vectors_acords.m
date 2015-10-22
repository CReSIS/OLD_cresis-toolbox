function create_vectors_acords(param,param_override)
% create_vectors_acords(param,param_override)
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
%  param = struct or .m file script which will create a struct
%  param_override = any fields here will be combined with param
%    using merge_structs
%
% Examples:
%   See default user settings inside this function.
%
% Authors: John Paden, Logan Smith
%
% See also: create_vectors_snow.m, plot_vectors.m,master.m,make_gps_*.m

% =====================================================================
% General Setup
% =====================================================================
%clear; % Useful when running as script
%close all; % Optional
tstart_create_vectors_acords = tic;
fprintf('\n\n==============================================\n\n');

% =====================================================================
% User Settings
% =====================================================================

if ~exist('param','var') || isempty(param)
  param = read_param_xls('E:\mcords_param_2010_Antarctica_DC8_P3.xls','20101013_04');
  
  % Input checking
  if ~exist('param','var')
    error('A struct of parameters must be passed in\n');
  end
  global gRadar;
  if exist('param_override','var')
    param_override = merge_structs(gRadar,param_override);
  else
    param_override = gRadar;
  end
end
param = merge_structs(param, param_override); 

% =====================================================================
% Load the GPS data
% =====================================================================

gps = load(ct_filename_support(param,param.vectors.gps.fn,'gps',true));

% The data only need 1 sample per second. This is a robust way
% to determine the sampling rate and to decimate by that.
decim_factor = max(1,round(1/median(diff(gps.gps_time))));
gps.gps_time = gps.gps_time(1:decim_factor:end);
gps.lat = gps.lat(1:decim_factor:end);
gps.lon = gps.lon(1:decim_factor:end);
gps.elev = gps.elev(1:decim_factor:end);

% Get the GPS seconds of day to sync to radar
utc_time_datenum = epoch_to_datenum(gps.gps_time-utc_leap_seconds(gps.gps_time(1)));
[year month day hour minute sec] = datevec(utc_time_datenum);

UTC_sod = (day-day(1))*86400+hour*3600+minute*60+sec;  % UTC seconds of day

% Remove repeat values in time
[UTC_sod sort_idxs] = unique(UTC_sod);
gps.lat = gps.lat(sort_idxs);
gps.lon = gps.lon(sort_idxs);
gps.elev = gps.elev(sort_idxs);
gps.gps_time = gps.gps_time(sort_idxs);
if(~isfield(gps,'time_offset'))
    gps.time_offset = 0;
end

% =====================================================================
% Get the files
% =====================================================================

if ~iscell(param.vectors.file.base_dir)
  param.vectors.file.base_dirs = {param.vectors.file.base_dir};
  param.vectors.file.adc_folder_name = {param.vectors.file.adc_folder_name};
  param.vectors.file.file_prefix = {param.vectors.file.file_prefix};
end

vectors_idx = 0;
for base_idx = 1:length(param.vectors.file.base_dirs)
  base_dir = param.vectors.file.base_dirs{base_idx};
  base_dir = ct_filename_data(param,base_dir);
  adc_folder_name = param.vectors.file.adc_folder_name{base_idx};
  if isfield(param.vectors.file,'file_prefix')
    file_prefix = param.vectors.file.file_prefix{base_idx};
  else
    file_prefix = '';
  end
  
  % Create sub-folder name for the particular receiver
  adc_idx_insert_idxs = strfind(adc_folder_name,'%d');
  mat_cmd = 'adc_folder_name = sprintf(adc_folder_name';
  for adc_idx_insert_idx = adc_idx_insert_idxs
    mat_cmd = [mat_cmd sprintf(', %d',param.vectors.file.adc)];
  end
  mat_cmd = [mat_cmd ');'];
  eval(mat_cmd);
  
  filepath = fullfile(base_dir, adc_folder_name);
  
  fprintf('Creating vectors for: %s\n', param.day_seg);
  file_mask = '';
  filenames = get_filenames(filepath, file_prefix, file_mask, '');
  if isempty(filenames)
    fprintf('No files match the mask:\n');
    fprintf('  path: %s\n', filepath);
    fprintf('  mask: *%s*.raw\n', file_mask);
    continue;
  end
    
  if param.vectors.file.stop_idx > length(filenames)
    stop_idx = length(filenames);
  else
    stop_idx = param.vectors.file.stop_idx;
  end
  file_idxs = param.vectors.file.start_idx:stop_idx;

  fprintf('  Input: %s\n', filepath);
  fprintf('  Output: %s\n', param.vectors.out_fn);
  
  for file_num = file_idxs
    fn = filenames{file_num};
    vectors_idx = vectors_idx + 1;

    % no need to determine version, ACORDS has one primary version/structure
    fname.name           = 'acords';
    fname.group         = '';
    fname.file_idx      = file_num;
    fname.radar_num     = 1;
    fname.rx            = 1;
    time_stamp          = fn(end-22:end-9);
    fname.datenum       = datenum(  str2double(time_stamp(1:4)), ...
                                    str2double(time_stamp(5:6)), str2double(time_stamp(7:8)), ...
                                    str2double(time_stamp(9:10)), str2double(time_stamp(11:12)), ...
                                    str2double(time_stamp(13:14)));
    fname.datevec       = datevec(fname.datenum);
    
    % Read first header
    [fid msg]       = fopen(fn,'r');
    if fid == -1
      fprintf('Unable to open file: %s\n',fn);
      error('File:Open',msg);
    end
    fseek(fid,4,'bof');
    
    seconds         = fread(fid,1,'int32');
    useconds        = fread(fid,1,'int32');
    comp_time       = seconds+(1e-6*useconds);
    
    hdr.time        = datenum(floor(datevec(comp_time/(24*60*60) + datenum(1970,1,1,0,0,0))));      % REPLACE WITH GIVEN SUBFUNCTION
    hdr.time        = hdr.time + gps.time_offset/(24*60*60);                                        % to account for computer time to gps time offset 
    [year month day hour minute sec] = datevec(hdr.time);
    hdr.sod(vectors_idx)         = hour*3600 + minute*60 + sec;                  % UTC seconds of day   
    fclose(fid);
  
    % Store the filename information
    fprintf('  File index %d (filename file_num %d) (%.1f sec)\n', ...
      file_num, fname.file_idx, toc(tstart_create_vectors_acords));
    vectors.file(vectors_idx).idx               = fname.file_idx;
    vectors.file(vectors_idx).rxs               = fname.rx;
    vectors.file(vectors_idx).recording_group   = fname.group;
    vectors.file(vectors_idx).base_dir          = base_dir;
    vectors.file(vectors_idx).adc_folder_name   = adc_folder_name;
    vectors.file(vectors_idx).datenum           = fname.datenum;
    vectors.fileNumber(vectors_idx)             = vectors_idx;
  end
end
  
% Check for seconds of day roll over and unwrap (assume jump backward of
% more than 23 hours is a roll over)
wrap_idxs = find(diff(hdr.sod) < -86300);
for wrap_idx = wrap_idxs
  hdr.sod(wrap_idx+1:end) = hdr.sod(wrap_idx+1:end) + 86400;
end
hdr.sod = hdr.sod + param.vectors.gps.time_offset;
         
% Synchronize to GPS data
vectors.lat = interp1(UTC_sod,gps.lat,hdr.sod);
vectors.lon = mod(interp1(UTC_sod,unwrap(gps.lon/180*pi),hdr.sod)*180/pi+180,360)-180;
vectors.elev = interp1(UTC_sod,gps.elev,hdr.sod);
vectors.gps_time = interp1(UTC_sod,gps.gps_time,hdr.sod);
vectors.gps_source = gps.gps_source;

% Make sure output directory exists
[out_dir out_name] = fileparts(ct_filename_support(param,param.vectors.out_fn,'vectors'));
if ~exist(out_dir,'dir')
  mkdir(out_dir);
end

param_vectors = param.vectors;
save(ct_filename_support(param,param.vectors.out_fn,'vectors'),'vectors','param_vectors');

return;
