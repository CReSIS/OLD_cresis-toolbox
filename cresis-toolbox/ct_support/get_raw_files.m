function [load_info,gps_time,recs] = get_raw_files(param,frm_id,imgs,rec_range,rec_range_type,out_dir)
% [load_info,gps_time,recs]= get_raw_files(param,frm_id,imgs,rec_range,rec_range_type)
%
% param: Can be either a string with the parameter spreadsheet filename OR
%  a struct containing radar, season, and optionally day_seg information.
%  .radar_name: string containing radar name (e.g. 'kuband2')
%  .season_name: string containing season name (e.g. '2012_Greenland_P3')
%  .day_seg: string containing the day segment (e.g. '20170412_01'). This
%    is not required if frm_id is a string with the day_seg in it.
% imgs: a cell array of wf-adc pair lists
% frm_id: One of these options:
%   1. string containing frame id (e.g. '20120514_01_317')
%   2. an integer containing the frame number (param.day_seg must be passed
%   in)
% rec_range: Specifies a range of records to load in. The units specified by
%   rec_range_type. Only the first and last element of rec_range are used.
% rec_range_type: String containing 'gps_time' or 'records'. Default is
%   'records'. If 'records' is used, either the param.day_seg field must be
%   defined or the frm_id must be a string with the day_seg in it.
% out_dir: Optional. May be left empty or not defined. Specifies an
%   output directory to copy raw files to.
%
% load_info: struct with file information for the range of data specified
%  .filenames: cell array of cells which contain the raw data filenames
%  .file_idx: cell array of integers which specify an index into .filenames
%    for each record
%  .offset: raw file byte offset to the beginning of each record
% gps_time: prints out the start and stop GPS times and puts them here
% recs: raw data records into csarp_support records file
%
% Example:
%   % Get filename information for a range of records
%   [load_info,gps_time,recs] = get_raw_files(struct('radar_name','rds','season_name','2017_Antarctica_TObas'),'20170122_01',{},[25092 27499])
%
%   % Get filename information for a particular frame
%   [load_info,gps_time,recs] = get_raw_files(struct('radar_name','accum2','season_name','2012_Greenland_P3'),'20120514_02_010');
%   [load_info,gps_time,recs] = get_raw_files('kuband_param_2012_Greenland_P3.xls','20120514_01_317');
%   [load_info,gps_time,recs] = get_raw_files(struct('radar_name','mcords2','season_name','2012_Greenland_P3'),'20120514_02_018');
%   load_info = get_raw_files(struct('radar_name','mcords','season_name','2011_Antarctica_DC8'),'20111030_02_002')
%   [load_info,gps_time,recs] = get_raw_files('accum_param_2017_Greenland_P3.xls','20170412_01_023',[],1.4920018728e9,'gps_time');
%
%   % Example copying files
%   load_info = get_raw_files(struct('radar_name','mcrds','season_name','2008_Greenland_TO'),'20080627_06_001','/tmp/)
%
% Author: John Paden
%
% See also: get_frame_id, get_raw_files.m, get_segment_file_list.m,
%   run_get_segment_file_list.m

param_fn = '';
if ischar(param)
  param_fn = ct_filename_param(param);
  clear param;
  if ~ischar(frm_id)
    error('param as a filename requires frm_id to be a frame ID string so the segment can be determined.');
  end
end
if ischar(frm_id)
  param.day_seg = frm_id(1:11);
elseif ~isfield(param,'day_seg')
  error('param.day_seg or frm_id as a frame id string must be provided');
end
if ischar(frm_id)
  frm = str2double(frm_id(end-2:end));
else
  frm = frm_id;
  frm_id = sprintf('%s_%03d', param.day_seg, frm);
end

if ~isempty(param_fn)
  % A filename was passed in for the radar parameters
  param = read_param_xls(ct_filename_param(param_fn),frm_id(1:11));
elseif isstruct(param)
  % A structure was passed in for the radar parameters
  [output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);
  param_fn = sprintf('%s_param_%s.xls', output_dir, param.season_name);
  param = read_param_xls(ct_filename_param(param_fn),frm_id(1:11));
else
  error('param must be a filename or structure.');
end

global gRadar;
param = merge_structs(gRadar,param);

% Populate imgs cell array with all wf-adc pairs if not specified
if isempty(imgs)
  for wf = 1:length(param.radar.wfs)
    for adc = 1:length(param.radar.wfs(wf).rx_paths)
      imgs{wf}(adc,1:2) = [wf adc];
    end
  end
end
% Populate wf_adc_list
wf_adc_list = [];
for img = 1:length(imgs)
  wf_adc_list(end+(1:size(imgs{img},1)),1:2) = imgs{img};
end

% Load the records file
records = records_load(param);

% Load the frames file
frames = frames_load(param);

%% Get the records associated with the frm or record
if exist('rec_range','var') && ~isempty(rec_range)
  if exist('rec_range_type','var') && strcmp(rec_range_type,'gps_time')
    % Use the provided gps_time range to determine the records to get
    start_rec = find(records.gps_time>=rec_range(1),1,'first');
    if isempty(start_rec)
      error('rec_range start time is after end of segment.');
    end
    stop_rec = find(records.gps_time>=rec_range(end),1,'first');
    if isempty(stop_rec)
      warning('rec_range stop time is after end of segment, setting to end of segment.');
      stop_rec = length(records.gps_time);
    end
  else
    % Use the provided record rnage to determine the records to get
    start_rec = rec_range(1);
    stop_rec = rec_range(end);
    good_recs = intersect(1:length(records.gps_time),start_rec:stop_rec);
    if numel(good_recs) < stop_rec-start_rec+1
      warning('Record range requests records that do not exist, truncating to valid records.');
      start_rec = good_recs(1);
      stop_rec = good_recs(end);
    end
  end

else
  % Use the provided frame ID to determine the records to get
  
  if frm > length(frames.frame_idxs)
    error('Frame %d > %d does not exist\n', frm, length(frames.frame_idxs));
  elseif frm == length(frames.frame_idxs)
    % Special case that handles last frame in segment
    start_rec = frames.frame_idxs(frm);
    stop_rec = length(records.lat);
  else
    start_rec = frames.frame_idxs(frm);
    stop_rec = frames.frame_idxs(frm+1)-1;
  end
end
recs = start_rec:stop_rec;

%% Get GPS Times
gps_time = records.gps_time(start_rec:stop_rec);

%% Get filenames
load_info = get_raw_files_sub(param,wf_adc_list,records,recs);

%% Reduce filenames list to just the files in the range specified
copy_fns = {};
for idx = 1:length(load_info.filenames)
  %% Copy files if out_dir specified
  if exist('out_dir','var') && ~isempty(out_dir)
    [base_dir,board_folder_name] = get_segment_file_list(param,load_info.board_idx(idx));
    for file_idx = 1:length(load_info.filenames{idx})
      copy_fn = load_info.filenames{idx}{file_idx};
      [~,copy_fn_name,copy_fn_ext] = fileparts(copy_fn);
      if ~any(strcmpi(copy_fn,copy_fns))
        % Some radars like MCRDS may list the same file multiple times
        out_fn = fullfile(out_dir,board_folder_name,[copy_fn_name copy_fn_ext]);
        out_fn_dir = fileparts(out_fn);
        fprintf('Copying %s\n  to %s (%s)\n', copy_fn, out_fn, datestr(now));
        if ~exist(out_fn_dir,'dir')
          mkdir(out_fn_dir);
        end
        copyfile(copy_fn,out_fn);
        copy_fns{end+1} = copy_fn;
      end
    end
  end
end

