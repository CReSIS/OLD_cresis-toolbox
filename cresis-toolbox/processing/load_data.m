function [data,hdr] = load_data(param)
% [data,hdr] = load_data(param)
%
% Function for loading data using records structure.
% The following options are available:
%   raw, pulse compressed, fast time decimated, burst EMI removal,
%   rxs, wfs, receiver combine, presumming
%
% param = struct controlling the loading and processing
%
% .load_data = structure controlling load_data processing
%  .radar_name = name of radar string
%  .season_name = name of season string
%  .day_seg = day-segment string
%  .records_fn = filename of records file
%  .recs = records to load (one indexed)
%  .imgs = cell vector of images to load, each image is Nx2 array of
%     wf/adc pairs
%     NOTE: wfs/adc pairs are not indices into anything, they are the absolute
%     waveform/adc numbers. The records file will be loaded to decode
%     which index each wf/adc belongs to.
%  .debug_level = debug level (scalar integer)
%  
%  load_data fields used by load_mcords_data.m (see that function for details)
%  .pulse_comp
%  .ft_dec
%  .ft_wind
%  .ft_wind_time
%  .presums
%  .combine_rx
%  .trim_vals
%  .pulse_rfi.en
%  .pulse_rfi.inc_ave
%  .pulse_rfi.thresh_scale
%  .radar
%   .Vpp_scale = radar Vpp for full scale quanitization
%   .rxs = struct array of receiver equalization coefficients
%    .wfs = struct array for each waveform
%     .chan_equal = scalar complex double (data DIVIDED by this)
%     .td = time delay correction (applied during pulse compression)
%  
% Example: run_load_data.m
%
% Author: John Paden
%
% See also run_load_data.m, load_mcords_data.m, pulse_compress.m

% =====================================================================
% Setup processing
% =====================================================================

try; toc; catch; tic; end;

physical_constants;

% Load records file
records_fn = ct_filename_support(param,param.load_data.records_fn,'records');
records = load(records_fn);
old_param_records = records.param_records;
recs(1) = param.load_data.recs(1);
recs(2) = param.load_data.recs(2);

% Create waveform parameters
if strcmpi(param.radar_name,'mcrds')
  [wfs,rec_data_size] = load_mcrds_wfs(records.settings, param, ...
    1:max(old_param_records.file.adcs), param.load_data);
  load_param.load.rec_data_size = rec_data_size;
elseif any(strcmpi(param.radar_name,{'acords','hfrds','mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
  [wfs,rec_data_size] = load_mcords_wfs(records.settings, param, ...
    1:max(old_param_records.records.file.adcs), param.load_data);
  load_param.load.rec_data_size = rec_data_size;
elseif strcmpi(param.radar_name,'icards')%try to add icards radar here--------------QISHI
  [wfs,rec_data_size] = load_icards_wfs(records.settings, param, ...
    1:max(old_param_records.records.file.adcs), param.load_data);
  load_param.load.rec_data_size = rec_data_size;
elseif any(strcmpi(param.radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3'}))
  [path name] = fileparts(records_fn);
  cdf_fn = fullfile(path, sprintf('%s.nc', name));
  try
    records.settings.nyquist_zone = ncread(cdf_fn,'settings(1).nyquist_zone',[1 load_param.load.recs(1)],[1 1]);
  end
  try
    records.settings.loopback_mode = ncread(cdf_fn,'settings(1).loopback_mode',[1 load_param.load.recs(1)],[1 1]);
  end
  wfs_idx = find(records.settings.wfs_records <= load_param.load.recs(1),1,'last');
  records.settings.wfs = records.settings.wfs(wfs_idx).wfs;
  wfs = load_fmcw_wfs(records.settings, param, ...
    1:max(old_param_records.records.file.adcs), param.load_data);
end

% Setup inputs for load_mcords_data
if ~isfield(param,'debug_level');
  load_param.debug_level = 1;
else
  load_param.debug_level = param.debug_level;
end

load_param.proc = param.load_data;

param.load.imgs = param.load_data.imgs;
% load_param.load.imgs = param.load_data.imgs;

load_param.wfs = wfs;

load_param.radar = param.radar;

load_param.load.rec_data_size = rec_data_size;

% =====================================================================
% Load data
% =====================================================================

% Determine where breaks in processing are going to occur
REC_BLOCK_SIZE = 10000;
if length(recs) < 2*REC_BLOCK_SIZE
  breaks = 1;
else
  breaks = 1:REC_BLOCK_SIZE:length(recs)-REC_BLOCK_SIZE;
end

global g_data;
g_data = [];

% Begin loading data
clear data;
for break_idx = 1:length(breaks)
  % Determine the current records being processed
  if break_idx < length(breaks)
    cur_recs = [recs(breaks(break_idx)) recs(breaks(break_idx+1)-1)];
  else
    cur_recs = [recs(breaks(break_idx)) recs(end)];
  end
  fprintf('Loading records %d to %d (%.1f sec)\n', ...
    cur_recs(1), cur_recs(end), toc);
  
  % =====================================================================
  % Load file information for each receiver channel
  % =====================================================================
  
  % Create a list of unique adcs required by the imgs list
  param.load.adcs = [];
  for idx = 1:length(param.load.imgs)
    new_adcs = abs(param.load.imgs{idx}(:,2:2:end)).';
    param.load.adcs = unique(cat(2, new_adcs(:).', param.load.adcs));
  end
  
  load_param.load.recs = cur_recs(1):cur_recs(end);
  
  if any(strcmpi(param.radar_name,{'icards','hfrds','mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
    % adc_headers: the actual adc headers that were loaded
    if ~isfield(old_param_records.records.file,'adc_headers') || isempty(old_param_records.records.file.adc_headers)
      old_param_records.records.file.adc_headers = old_param_records.records.file.adcs;
    end
    
    % boards_headers: the boards that the actual adc headers were loaded from
    boards_headers = adc_to_board(param.radar_name,old_param_records.records.file.adcs);
    
    for idx = 1:length(param.load.adcs)
      % adc: the specific ADC we would like to load
      adc = param.load.adcs(idx);
      % adc_idx: the records file index for this adc
      adc_idx = find(old_param_records.records.file.adcs == adc);
      if isempty(adc_idx)
        error('ADC %d not present in records file\n', adc);
      end
      
      % board: the board associated with the ADC we would like to load
      board = adc_to_board(param.radar_name,adc);
      % board_header: the board headers that we will use with this ADC
      board_header = adc_to_board(param.radar_name,old_param_records.records.file.adc_headers(adc_idx));
      % board_idx: the index into the records board list to use
      board_idx = find(board_header == boards_headers);
      
      % Just get the file-information for the records we need
      load_param.load.file_idx{idx} = relative_rec_num_to_file_idx_vector( ...
        load_param.load.recs,records.relative_rec_num{board_idx});
      load_param.load.offset{idx} = records.offset(board_idx,load_param.load.recs);
      file_idxs = unique(load_param.load.file_idx{idx});
      
      % Recognize if first record is really from previous file and it is a
      % valid record (i.e. offset does not equal -2^31)
      if sign(load_param.load.offset{idx}(1)) < 0 && load_param.load.offset{idx}(1) ~= -2^31
        file_idxs = [file_idxs(1)-1 file_idxs];
      end
      
      % Just copy the filenames we need
      load_param.load.filenames{idx}(file_idxs) = records.relative_filename{board_idx}(file_idxs);

      % Modify filename according to channel
      for file_idx = 1:length(load_param.load.filenames{idx})
        if any(strcmpi(param.radar_name,{'mcords5'}))
          load_param.load.filenames{idx}{file_idx}(9:10) = sprintf('%02d',board);
        end
      end
      
      filepath = get_segment_file_list(param,adc);
      
      % Convert relative file paths into absolute file paths if required,
      % also corrects filesep (\ and /)
      for file_idx = 1:length(load_param.load.filenames{idx})
        load_param.load.filenames{idx}{file_idx} ...
          = fullfile(filepath,load_param.load.filenames{idx}{file_idx});
      end
    end
    load_param.load.file_version = param.records.file_version;
  elseif strcmpi(param.radar_name,'mcrds')
    load_param.load.file_rec_offset = records.file_rec_offset;
    load_param.load.filenames = records.filenames;
    base_dir = ct_filename_data(ct_filename_param,param.get_heights.file.base_dir);
    adc_folder_name = param.get_heights.file.adc_folder_name;
    load_param.load.filepath = fullfile(base_dir, adc_folder_name);
    load_param.load.wfs = records.wfs;
    load_param.load.wfs_file = records.wfs_file;
  elseif strcmpi(param.radar_name,'acords')
    load_param.load.file_rec_offset = records.relative_rec_num;
    load_param.load.filenames = records.relative_filename;
    base_dir = ct_filename_data(param,param.vectors.file.base_dir);
    adc_folder_name = param.vectors.file.adc_folder_name;
    load_param.load.filepath = fullfile(base_dir, adc_folder_name);
    load_param.load.wfs = records.settings.wfs;
    load_param.load.file_version = param.records.file_version;
    load_param.load.offset = records.offset;
    load_param.load.wfs_records = records.settings.wfs_records;

%   elseif strcmpi(param.radar_name,'icards')%try to add a situation of icards----------QISHI
%       load_param.load.file_rec_offset = records.relative_rec_num;%file_rec_offset---->relative_rec_num
%       load_param.load.filenames = records.relative_filename;     %filenames------>relative_filename
%       % %     base_dir = ct_filename_data(ct_filename_param,param.vectors.file.base_dir);%get_heights---->vectors
%       base_dir = param.vectors.file.base_dir;
%       adc_folder_name = param.vectors.file.adc_folder_name;                      %get_heights---->vectors
%       load_param.load.filepath = fullfile(base_dir, adc_folder_name);
%       load_param.load.wfs = records.settings.wfs;          %records.wfs---->records.settings.wfs
%       load_param.load.wfs_file = records.settings.wfs_file;%records.wfs_file---->records.settings.wfs_file(wfs_file is always assumed to be 1 for icards)

  elseif any(strcmpi(param.radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3'}))
    load_param.load.offset = records.offset;
    load_param.load.file_rec_offset = records.relative_rec_num;
    load_param.load.filenames = records.relative_filename;
    base_dir = ct_filename_data(param,param.vectors.file.base_dir);
    adc_folder_name = param.vectors.file.adc_folder_name;
    load_param.load.filepath = fullfile(base_dir, adc_folder_name);
    load_param.load.wfs = records.settings.wfs;
    load_param.load.radar_name = param.radar_name;
    load_param.load.season_name = param.season_name;
    load_param.load.tmp_path = param.tmp_path;
    load_param.load.file_version = param.records.file_version;
    %   if any(strcmpi(param.radar_name,{'snow','kuband'}))
    %     load_param.load.header_size = 4*40;
    %   elseif any(strcmpi(param.radar_name,{'snow2','kuband2'}))
    %     load_param.load.header_size = 32;
    %   elseif any(strcmpi(param.radar_name,{'snow3','kuband3'}))
    %     if param.records.file_version == 4
    %       load_param.load.header_size = 32;
    %     else
    %       load_param.load.header_size = 48;
    %     end
    %   end
  end
  
  % =====================================================================
  % Setup control parameters for load_mcords_data
  
  load_param.load.imgs = param.load.imgs;
  load_param.load.adcs = param.load.adcs;

  % =====================================================================
  % Load data
  % =====================================================================
  
  % Load data into g_data
  if strcmpi(param.radar_name,'mcords')
    load_mcords_data(load_param);
  elseif any(strcmpi(param.radar_name,{'hfrds','mcords2','mcords3','mcords4','mcords5'}))
    load_mcords2_data(load_param);
  elseif strcmpi(param.radar_name,'mcrds')
    load_mcrds_data(load_param);
  elseif strcmpi(param.radar_name,'acords')
    load_acords_data(load_param);
  elseif any(strcmpi(param.radar_name,{'icards'}))
    load_icards_data(load_param);%try to add a loader of icards----------QISHI
  elseif strcmpi(param.radar_name,'accum2')
    load_accum2_data(load_param);
  elseif any(strcmpi(param.radar_name,{'kuband','snow','kuband2','snow2','kuband3','snow3'}))
    wfs.time = load_fmcw_data(load_param);
  end
  
  % =====================================================================
  % Concatenate data
  % =====================================================================
  if break_idx == 1
    data = g_data;
  else
    for img_idx = 1:length(g_data)
      data{img_idx} = cat(2, data{img_idx}, g_data{img_idx});
    end
  end
  
end

% =====================================================================
% Create output hdr struct
% =====================================================================
hdr.records = fir_dec(recs(1):recs(2),param.load_data.presums);
hdr.param = param;
if isfield(records,'lat')
  hdr.lat = fir_dec(records.lat(recs(1):recs(2)),param.load_data.presums);
  hdr.lon = fir_dec(records.lon(recs(1):recs(2)),param.load_data.presums);
  hdr.elev = fir_dec(records.elev(recs(1):recs(2)),param.load_data.presums);
  hdr.roll = fir_dec(records.roll(recs(1):recs(2)),param.load_data.presums);
  hdr.pitch = fir_dec(records.pitch(recs(1):recs(2)),param.load_data.presums);
  hdr.heading = fir_dec(records.heading(recs(1):recs(2)),param.load_data.presums);
  hdr.gps_time = fir_dec(records.gps_time(recs(1):recs(2)),param.load_data.presums);
  hdr.gps_source = records.gps_source;
end
if isfield(records,'surface')
    hdr.surface = fir_dec(records.surface(recs(1):recs(2)),param.load_data.presums);
end
hdr.wfs = wfs;

return;
