function create_vectors_mcrds(param,param_override)
% create_vectors_mcrds(param,param_override)
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
% Authors: John Paden
%
% See also: create_vectors_snow.m, plot_vectors.m,master.m,make_gps_*.m

% =====================================================================
% General Setup
% =====================================================================

dbstack_info = dbstack;
if ~exist('param','var') || isempty(param) || length(dbstack_info) == 1
  % =====================================================================
  % Debug Setup
  % =====================================================================
  param = read_param_xls(ct_filename_param('rds_param_2006_Greenland_TO.xls'),'20060526_01');
  
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
% Get the list of files to include in this vectors file
% =====================================================================
[base_dir,~,fns,file_idxs] = get_segment_file_list(param);

clear vectors;
vectors_idx = 0;
for file_num = 1:length(file_idxs)
  fn = fns{file_idxs(file_num)};
  
  vectors_idx = vectors_idx + 1;
  
  % Store the file number index
  vectors.fileNumber(vectors_idx) = file_num;
  
  % Build up structure in fname to append to file struct vector
  [fn_dir fn_name] = fileparts(fn);
  fname.file_prefix = fn_name(1:find(fn_name=='.',1));
  fname.idx = str2double(fn(end-7:end-4));
  time_stamp          = fn(end-22:end-9);
  fname.datenum       = datenum(  str2double(time_stamp(1:4)), ...
    str2double(time_stamp(5:6)), str2double(time_stamp(7:8)), ...
    str2double(time_stamp(9:10)), str2double(time_stamp(11:12)), ...
    str2double(time_stamp(13:14)));
  fname.base_dir          = base_dir;
  fname.adc_folder_name   = param.vectors.file.adc_folder_name;
  vectors.file(vectors_idx) = fname;
  
  % Read first header
  Header = basic_load_mcrds_hdr(fn);
  
  % Read in the computer and radar time for the first record in the file
  [fid msg]       = fopen(fn,'r');
  if fid == -1
    fprintf('Unable to open file: %s\n',fn);
    error('File:Open',msg);
  end
  % Skip header and first DAQ bit set
  num_of_header_bytes = 40 + 444 + 176*Header.NumberWaveforms;
  status          = fseek(fid,num_of_header_bytes+4,'bof');
  % Read in 32 bit seconds since Jan 1, 1970 and 32 bit microseconds field
  vectors.comp_time(vectors_idx) = fread(fid,1,'int32=>double') + fread(fid,1,'int32=>double')*1e-6;
  % Read in 10 MHz counts field
  vectors.radar_time(vectors_idx) = fread(fid,1,'int64=>double') * 1e-7;
  fclose(fid);
  
  % Print basic file information to console
  fprintf('  File index %d (filename file_num %d) (%s sec)\n', ...
    file_num, fname.idx, datestr(now));
end

if param.vectors.gps.time_offset ~= 0
  warning('MCRDS usually has 0 time offset');
end
vectors.radar_time = vectors.radar_time + param.vectors.gps.time_offset;

%% Synchronize to GPS data
vectors = sync_radar_to_gps(param,vectors,vectors.radar_time,vectors.comp_time);

%% Make sure output directory exists
[out_dir out_name] = fileparts(ct_filename_support(param,param.vectors.out_fn,'vectors'));
if ~exist(out_dir,'dir')
  mkdir(out_dir);
end

%% Save outputs
param_vectors = param.vectors;
save(ct_filename_support(param,param.vectors.out_fn,'vectors'),'vectors','param_vectors');

return;
