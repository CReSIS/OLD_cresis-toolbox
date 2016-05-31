function create_records_icards(param,param_override)
%CHANGE ALL DESCRIPTION
% create_records_icards(param,param_override)
%
% Corrects jumps in utc time, typically one second jumps. This script
% obtains headers from data drives indicated in the param
% spreadsheets. After loading the files using basic_load.m the header files
% are saved in the support directories for the specific radar in .mat form
% for quicker access.
%
% Can be run as a function by passing in the param argument
% or a script (by setting the default value of param).
%
% Output file contains:
% hdr: structure with the following fields
%    .utc_time_sod: vector of corrected utc times
%    .seconds: vector
%    .fraction: vector
%
% Inputs:
% param = struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Authors: Aric Beaver, John Paden
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
  param = read_param_xls(ct_filename_param('rds_param_2002_Greenland_P3.xls'),'20020520_01');
  param.vectors.file.base_dir='Z:\ICARDS\2002\';  
  
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
%% Check input parameters
% =====================================================================

if ~isfield(param.records,'force_all') || isempty(param.records.force_all)
  param.records.force_all = false;
end

if ~isfield(param.records,'file_version') || isempty(param.records.file_version)
  param.records.file_version = 1;
end

if ~isfield(param.records,'debug_level') || isempty(param.records.debug_level)
  param.records.debug_level = 1;
end

% =====================================================================
%% Setup the scheduler
% =====================================================================

param.sched.type = 'no scheduler';
% =====================================================================
%% Get the list of raw files to load in
% =====================================================================
% =====================================================================
%% Load headers from radar data files
% =====================================================================
param.file_regexp = '\S+\.[0-9][0-9][0-9]$';
full_dir = fullfile(param.vectors.file.base_dir,param.vectors.file.adc_folder_name);
fns = get_filenames(full_dir,'','','',struct('regexp',param.file_regexp));
valid_data_file=icards_data_ignore_list(fns,full_dir);
fns=fns(valid_data_file);
file_idxs=param.vectors.file.start_idx:param.vectors.file.stop_idx;
param_1=param;
param_1.day_seg=[];
hdrs.offset=[];

if file_idxs(1)==1
  first_seg=1;
else
  first_seg=0;
  file_idx_previous=file_idxs(1)-1;
end


num_rec_sample=0;%initial number of samples per record
hdrs.nmea_time=[];
for file_idxs_idx = 1:length(file_idxs)
  file_idx = file_idxs(file_idxs_idx);
  fn=fns{file_idx};
  [~,fn_name,fn_ext] = fileparts(fn);
  fprintf('  File %s %d of %d (%s)\n', fn_name, file_idxs_idx, length(file_idxs), datestr(now));
  tmp_hdr_fn = ct_filename_tmp(param_1,'','headers',[fn_name fn_ext '.mat']);
  tmp_hdr_fn_dir = fileparts(tmp_hdr_fn);
  if ~exist(tmp_hdr_fn_dir,'dir')
    mkdir(tmp_hdr_fn_dir);
  end
  hdr = load(tmp_hdr_fn); 
  
  if (all(isnan(hdr.nmea_time)))||(all(isnan(hdr.nmea_lat)))||(all(isnan(hdr.nmea_lon)))||(all(isnan(hdr.nmea_elev)))%to ingnore this temporary file which contains only NaN---qishi
    warning('no valid data in %s except NaN,to load next file if there are any\n',fn);
    continue;
  end
  
  if ~first_seg%if this is not the first segment,we need load the last time point of previous segment for interpolation to avoid "out of records" warning----qishi
    fn_previous=fns{file_idx_previous};
    [~,fn_name_previous,fn_ext_previous] = fileparts(fn_previous);
    tmp_hdr_fn_previous = ct_filename_tmp(param_1,'','headers',[fn_name_previous fn_ext_previous '.mat']);
    hdr_previous=load(tmp_hdr_fn_previous);
    while (all(isnan(hdr_previous.nmea_time)))%we need to find the last time point previous file when the file does not contain NaN only---qishi
      file_idx_previous=file_idx_previous-1;
      [~,fn_name_previous,fn_ext_previous] = fileparts(fn_previous);
      tmp_hdr_fn_previous = ct_filename_tmp(param_1,'','headers',[fn_name_previous fn_ext_previous '.mat']);
      hdr_previous=load(tmp_hdr_fn_previous);
    end
  
    if hdr.nmea_time(1)<hdr_previous.nmea_time(end)
      warning('the first nmea time of this file is smaller than the last time point of prvious file, correct this in interpolation step\n');
      previous_mark=1;
      previous_time=hdr_previous.nmea_time(end);
    else
      previous_mark=0;
      previous_time=[];
    end
  else
    previous_mark=0;
    previous_time=[];
  end
  
  
  if isempty(hdrs.nmea_time)
    fisrt_valid_mark=0;
  else
    fisrt_valid_mark=1;
  end
  
  if (file_idxs_idx == 1)||(fisrt_valid_mark==0)
    hdrs.filenames{1}=[fn_name fn_ext];
    hdrs.file_rec_offset = 1;
    hdrs.nmea_time = hdr.nmea_time;
  else
    hdrs.filenames{file_idxs_idx}=[fn_name fn_ext];
    hdrs.file_rec_offset(file_idxs_idx) = length(hdrs.nmea_time) + 1;
    hdrs.nmea_time = cat(2,hdrs.nmea_time,hdr.nmea_time);
  end
%==========================================================================  
  %create offset parameter here: offset vector contains the location of
  %first I sample from every record. These hearder information is observed
  %from function "icards_get_data" with "Z:\ICARDS\2002\may18\may18_02.074"
  %as an example. The data structure may be different for other day's data.
  %so this script should be modified!!??-------qishi/20150825
  num_records=length(hdr.nmea_time);
  header_size=64;
  secondary_header_size=12;
  if num_rec_sample==0
     num_rec_sample=size(icards_get_data(fn,2),1);  %if>0,coherent data in this segment
     if num_rec_sample==0
       num_rec_sample=size(icards_get_data(fn,1),1);%incoherent data in this segment
       fprintf('Incoherent data stored in this file \n');
     else
       fprintf('Coherent data stored in this file \n');
     end
  else
     num_rec_sample=num_rec_sample;
  end
  sample_size=2;
  rec_data_size=num_rec_sample*sample_size;
  hdrs.offset=[hdrs.offset header_size+secondary_header_size+rec_data_size*(0:num_records-1)];%this offset is the location of each I sample
%==========================================================================  
end
hdrs.filenames=hdrs.filenames(logical(hdrs.file_rec_offset));% arrange filenames
hdrs.file_rec_offset=hdrs.file_rec_offset(logical(hdrs.file_rec_offset));% arrange rec offset
% =====================================================================
%% interpolation to make time monotonically increasing 
% =====================================================================
if any(isnan(hdrs.nmea_time))%deal with possible NaNs in time sequence
  warning('NaN found in this file while creating records of this segment,correcting now\n');
  hdrs.nmea_time=create_records_icards_dealwithNaN(hdrs.nmea_time);
end
hdrs.nmea_time=create_records_icards_interpolation(hdrs.nmea_time,previous_mark,previous_time);% fix time sequence---qishi
hdr.nmea_time=hdrs.nmea_time;
hdr.filenames=hdrs.filenames;
hdr.file_rec_offset = hdrs.file_rec_offset;
hdr.offset=hdrs.offset;
%% Save workspace in case there is a failure
if isfield(param,'tmp_path') && ~isempty(param.tmp_path)
  fn = ct_filename_tmp(param,param.records.records_fn,'records','workspace');
  fprintf('Saving workspace %s (%s)\n', fn, datestr(now));
  try
    fn_dir = fileparts(fn);
    if ~exist(fn_dir,'dir')
      mkdir(fn_dir);
    end
    save(fn);
  catch
    fprintf('  Saving workspace failed\n');
  end
end
 create_records_icards_sync;

return;

