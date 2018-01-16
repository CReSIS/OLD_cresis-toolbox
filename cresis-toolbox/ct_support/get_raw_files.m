function [load_info,gps_time,recs] = get_raw_files(param,frm_id,out_dir,rec_range,rec_range_type)
% [load_info,gps_time,recs]= get_raw_files(param,frm_id,out_dir,rec_range,rec_range_type)
%
% param: Can be either a string with the parameter spreadsheet filename OR
%  a struct containing radar, season, and optionally day_seg information.
%  .radar_name: string containing radar name (e.g. 'kuband2')
%  .season_name: string containing season name (e.g. '2012_Greenland_P3')
%  .day_seg: string containing the day segment (e.g. '20170412_01'). This
%    is not required if frm_id is a string with the day_seg in it.
% frm_id: One of these options:
%   1. string containing frame id (e.g. '20120514_01_317')
%   2. an integer containing the frame number (param.day_seg must be passed
%   in)
% out_dir: Optional. May be left empty or not defined. Specifies an
%   output directory to copy files to.
% rec_range: Specifies a range of records to load in the used specified by
%   rec_range_type. Only the first and last element of rec_range are used.
%   If this method is used, either the param.day_seg field must be defined
%   or the frm_id must be a string with the day_seg in it.
% rec_range_type: String containing 'gps_time' or 'records'. Default is
%   'records'.
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
%   fns = get_raw_files(struct('radar_name','mcords3','season_name','2014_Greenland_P3'),'20140401_03',[],[10250 10250])
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

if ischar(param)
  param_fn = param;
  param_fn = ct_filename_param(param);
  if ischar(frm_id)
    day_seg = frm_id(1:11);
  else
    error('param as a filename requires frm_id to be a frame ID string');
  end
elseif ischar(frm_id)
  param.day_seg = frm_id(1:11);
elseif ~isfield(param,'day_seg')
  error('param.day_seg or frm_id as a frame id string must be provided');
end
if ischar(frm_id)
  frm = str2double(frm_id(end-2:end));
else
  frm = frm_id;
end

if ischar(param)
  param = read_param_xls(ct_filename_param(param_fn),frm_id(1:11));
else
  [output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);
  param_fn = sprintf('%s_param_%s.xls', output_dir, param.season_name);
  param = read_param_xls(ct_filename_param(param_fn),frm_id(1:11));
end

global gRadar;
param = merge_structs(gRadar,param);

% Load the records file
records_fn = ct_filename_support(param,'','records');
% load(records_fn);
records = load(records_fn);

% Load the frames file
frames_fn = ct_filename_support(param,'','frames');
load(frames_fn);

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
records.offset = records.offset(:,recs);
load_info = get_raw_files_sub(param,records,recs);

%% Reduce filenames list to just the files in the range specified
copy_fns = {};
for chan = 1:length(load_info.filenames)
  [file_idxs,idx_mapping1,idx_mapping2] = unique(load_info.file_idx{chan});
  load_info.filenames{chan} = load_info.filenames{chan}(file_idxs);
  load_info.file_idx{chan} = idx_mapping2;
  %% Copy files if out_dir specified
  if exist('out_dir','var') && ~isempty(out_dir)
    [base_dir,adc_folder_name,~,~] ...
      = get_segment_file_list(param,param.records.file.adcs(chan),true);
    for file_idx = 1:length(load_info.filenames{chan})
      copy_fn = load_info.filenames{chan}{file_idx};
      [~,copy_fn_name,copy_fn_ext] = fileparts(copy_fn);
      if ~any(strcmpi(copy_fn,copy_fns))
        % Some radars like MCRDS may list the same file multiple times
        out_fn = fullfile(out_dir,adc_folder_name,[copy_fn_name copy_fn_ext]);
        out_fn_dir = fileparts(out_fn);
        fprintf('Copying %s\n  to %s\n', copy_fn, out_fn);
        if ~exist(out_fn_dir,'dir')
          mkdir(out_fn_dir);
        end
        copyfile(copy_fn,out_fn);
        copy_fns{end+1} = copy_fn;
      end
    end
  end
end

