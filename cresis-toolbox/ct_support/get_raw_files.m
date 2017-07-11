function [fns,gps_time] = get_raw_files(param,frm_id,out_dir)
% [fns,gps_time]= get_raw_files(param,frm_id,out_dir)
%
% param = struct containing radar and season information
%  .radar_name = string containing radar name (e.g. 'kuband2')
%  .season_name = string containing season name (e.g. '2012_Greenland_P3')
% frm_id = string containing frame id (e.g. '20120514_01_317')
%
% fns = cell structure with strings containing the raw filenames
%       associated with this frame
%   For mcords2:
%     fns is a 1xB cell vector of Nx1 cell vectors of strings where B is
%     the number of boards and N is the number of files from each board
%   For everything else:
%     fns is an Nx1 cell vector of strings
% gps_time = prints out the start and stop GPS times and puts them here
%
% Example:
%   [fns,gps_time] = get_raw_files(struct('radar_name','kuband2','season_name','2012_Greenland_P3'),'20120514_01_317');
%   fns = get_raw_files(struct('radar_name','accum2','season_name','2012_Greenland_P3'),'20120514_02_010')
%   fns = get_raw_files(struct('radar_name','mcords2','season_name','2012_Greenland_P3'),'20120514_02_018')
%   fns = get_raw_files(struct('radar_name','mcords','season_name','2011_Antarctica_DC8'),'20111030_02_002')
%   fns = get_raw_files(struct('radar_name','mcrds','season_name','2008_Greenland_TO'),'20080627_03_001')
%
%   fns = get_raw_files('C:\Users\dangermo\Documents\scripts\branch\params-cr1\rds_param_2011_Greenland_P3.xls','20110506_01_018','K:\paden')
%
% Author: John Paden

if ischar(param)
  param = read_param_xls(param,frm_id(1:11));
end

global gRadar;
param = merge_structs(gRadar,param);
param.day_seg = frm_id(1:11);

[output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);

% Load the records file
records_fn = ct_filename_support(param,'','records');
% load(records_fn);
records = load(records_fn);

% Load the frames file
frames_fn = ct_filename_support(param,'','frames');
load(frames_fn);

frm = str2double(frm_id(end-2:end));

%% Get the records associated with frm
if frm > length(frames.frame_idxs)
  error('Frame does not exist\n');
elseif frm == length(frames.frame_idxs)
  % Special case that handles last frame in segment
  start_rec = frames.frame_idxs(frm);
  stop_rec = length(records.lat);
else
  start_rec = frames.frame_idxs(frm);
  stop_rec = frames.frame_idxs(frm+1)-1;
end

%% Get filenames associated with the records in the frame
if any(strcmpi(radar_name,{'kuband','kuband2','kuband3','snow','snow2','snow3','snow5','snow8'}))
  % Get the raw file numbers associated with the frm
  start_filenum = find(records.relative_rec_num{1} <= start_rec,1,'last');
  stop_filenum = find(records.relative_rec_num{1} <= stop_rec,1,'last');
  % Get the filenames associated with the frame
  fns = records.relative_filename{1}(start_filenum:stop_filenum).';
elseif any(strcmpi(radar_name,{'accum2'}))
  % Get the raw file numbers associated with the frm
  start_filenum = find(records.relative_rec_num{1} <= start_rec,1,'last');
  stop_filenum = find(records.relative_rec_num{1} <= stop_rec,1,'last');
  % Get the filenames associated with the frame
  fns = records.relative_filename{1}(start_filenum:stop_filenum).';
elseif any(strcmpi(radar_name,{'acords','mcrds'}))
  % Get the raw file numbers associated with the frm
  start_filenum = find(records.relative_rec_num{1} <= start_rec,1,'last');
  stop_filenum = find(records.relative_rec_num{1} <= stop_rec,1,'last');
  % Get the filenames associated with the frame
  fns = records.relative_filename{1}(start_filenum:stop_filenum).';
elseif any(strcmpi(radar_name,{'mcords','mcords2','mcords3','mcords4','mcords5'}))
  % Get the filenames associated with the frame (channel 0, 1, 2, etc)
  filenames_idxs = 1:length(records.relative_filename);
  for fns_idx = 1:length(filenames_idxs)
    filenames_idx = filenames_idxs(fns_idx);
    % Get the raw file numbers associated with the frm
    %   start_filenum = records.relative_rec_num{filenames_idx}(start_rec);
    %   stop_filenum = records.relative_rec_num{filenames_idx}(stop_rec);
    start_filenum = find(records.relative_rec_num{filenames_idx} <= start_rec,1,'last');
    stop_filenum = find(records.relative_rec_num{filenames_idx} <= stop_rec,1,'last');
    fns{fns_idx} = records.relative_filename{filenames_idx}(start_filenum:stop_filenum);
  end
else
  error('Radar name %s not supported by this function yet', param.radar_name);
end

%% Get GPS Times
fprintf('GPS time from %s to %s\n', ...
  datestr(epoch_to_datenum(records.gps_time(start_rec))), ...
  datestr(epoch_to_datenum(records.gps_time(stop_rec))))
gps_time = records.gps_time(start_rec:stop_rec);

if exist('out_dir','var')
  if any(strcmpi(radar_name,{'mcords2','mcords3'}))
    board = [];
    for adc_idx = 1:length(param.records.file.adcs)
      adc = param.records.file.adcs(adc_idx);
      board(adc_idx) = adc_to_board(param.radar_name,adc);
    end
    [boards boards_idx] = unique(board);
    for board_idx = 1:length(boards)
      [base_dir,adc_folder_name,copy_fns,file_idxs] ...
        = get_segment_file_list(param,param.records.file.adcs(boards_idx(board_idx)));
      start_filenum = find(records.relative_rec_num{1} <= start_rec,1,'last');
      stop_filenum = find(records.relative_rec_num{1} <= stop_rec,1,'last');
      for file_idx = start_filenum:stop_filenum
        copy_fn = copy_fns{file_idx};
        [~,copy_fn_name,copy_fn_ext] = fileparts(copy_fn);
        out_fn = fullfile(out_dir,adc_folder_name,[copy_fn_name copy_fn_ext]);
        out_fn_dir = fileparts(out_fn);
        fprintf('Copying %s to\n  %s\n', copy_fn, out_fn);
        if ~exist(out_fn_dir,'dir')
          mkdir(out_fn_dir);
        end
        copyfile(copy_fn,out_fn);
      end
    end
  end
end

end
