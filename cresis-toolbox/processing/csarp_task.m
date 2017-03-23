function [success] = csarp_task(param)
% [success] = csarp_task(param)
%
% SAR process a chunk of data (a frame is divided into chunks
% to keep the memory requirements low enough).  This function
% is called from csarp.m.
%
% param = structure controlling the SAR processor
%  .debug_level = scalar integer, debugging level
%  .radar = structured used by load_mcords_hdr
%
%  .load = structure containing info about what data to load
%   .records_fn = records filename
%   .recs = 2 element vector containing the start and stop records
%   .imgs = cell vector of images to load, each image is Nx2 array of
%     wf/adc pairs
%     NOTE: wfs/ads are not indices into anything, they are the absolute
%     waveform/adc numbers. The records file will be loaded to decode
%     which index each wf/adc belongs to.
%
%  .proc = structure about which frame to process and how it is broken
%    into chunks
%   .frm = only used to determine the file name of the output
%   .output_along_track_offset = along-track offset from first input
%      record to the first output range line
%   .output_along_track_Nx = length of output in along-track
%
%  .csarp = structure containing SAR processing parameters
%   .file = struct with input file information
%      .base_dir: string, e.g. '/cresis/data3/MCoRDS/2011_Greenland_P3/'
%      .adc_folder_name = string, e.g. '20110507/board%b/seg_01'
%      .file_prefix = string, e.g. ''
%   .out_path = output path string
%   .combine_rx = boolean, combine channels before SAR processing
%   .coh_noise_removal = boolean, slow-time DC removal
%   .lever_arm_fh = string containing function name to lever arm
%   .mocomp = struct for motion compensation
%      .en = boolean, apply motion compensation
%      .type = see motion_comp.m for details
%      .uniform_en = boolean, resample data to uniform sampling in along-track
%   .sar_type = string, 'f-k' or 'tdc'
%   .sigma_x = along-track output sample spacing (meters)
%   .sub_aperture_steering = vector of doppler centroids to process to
%     (i.e. subapertures) normalized to the doppler bandwidth
%   .st_wind = function handle for slow time decimation
%   .start_eps = epsilon to use for sub-surface
%
% Fields used by load_mcords_data.m (see that function for details)
%  .pulse_rfi.en
%  .pulse_rfi.inc_ave
%  .pulse_rfi.thresh_scale
%  .trim_vals
%  .pulse_comp
%  .ft_dec
%  .ft_wind
%  .ft_wind_time
%  .radar.rx_path.chan_equal
%  .radar.rx_path.td
%
% success = boolean which is true when the function executes properly
%   if a task fails before it can return success, then success will be
%   empty
%
% Authors: William Blake, John Paden
%
% See also: csarp.m
%

%% Initialization and checking arguments
physical_constants;

global g_data;

[output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);

param.load.records_fn = ct_filename_support(param,param.records.records_fn,'records');

if param.csarp.combine_rx && param.csarp.mocomp.en
  warning('CSARP motion compensation mode must be 0 for combine_rx (setting to 0)');
  param.csarp.mocomp.en = 0;
end

if ~isfield(param.csarp,'start_eps') || isempty(param.csarp.start_eps) || param.csarp.start_eps == 0
  param.csarp.start_eps = er_ice;
end

if ~isfield(param.csarp,'trim_vals') || isempty(param.csarp.trim_vals)
  param.csarp.trim_vals = [0 0];
end

if ~isfield(param.csarp,'coh_noise_removal') || isempty(param.csarp.coh_noise_removal) ...
    || ~param.csarp.coh_noise_removal
  default_coh_noise_method = 0;
elseif param.csarp.coh_noise_removal
  default_coh_noise_method = 1;
end

if ~isfield(param.csarp,'coh_noise_method') || isempty(param.csarp.coh_noise_method)
  param.csarp.coh_noise_method = default_coh_noise_method;
end

if ~isfield(param.csarp,'coh_noise_arg')
  param.csarp.coh_noise_arg = [];
end

if ~isfield(param.csarp,'pulse_rfi') || isempty(param.csarp.pulse_rfi)
  param.csarp.pulse_rfi.en = 0;
end
 
if ~isfield(param.records,'file_version')
  param.records.file_version = [];
end

if ~isfield(param.csarp.mocomp,'filter')
  param.csarp.mocomp.filter = '';
end

if ~isfield(param.csarp,'ground_based')
  param.csarp.ground_based = 0;
end

if ~isfield(param.csarp,'presums')
  param.csarp.presums = 1;
end

if ~isfield(param.csarp,'deconvolution') || isempty(param.csarp.deconvolution)
  param.csarp.deconvolution = 0;
end

if ~isfield(param.csarp,'psd_smooth') || isempty(param.csarp.psd_smooth)
  param.csarp.psd_smooth = 0;
end

%% Load record information
% =====================================================================
load_param.load.recs = param.load.recs;
orig_records = read_records_aux_files(param.load.records_fn,param.load.recs);
old_param_records = orig_records.param_records;
old_param_records.gps_source = orig_records.gps_source;
start_time_for_fn = orig_records.gps_time(1);

%% Get the new surface
% =====================================================================
if isfield(param.csarp,'surface_src') && ~isempty(param.csarp.surface_src)
  
  % Get the generic layer data path
  layer_path = fullfile(ct_filename_out(param,param.csarp.surface_src,'',0));
  
  % Load the current frame
  layer_fn = fullfile(layer_path,sprintf('Data_%s_%03d.mat',param.day_seg,param.proc.frm));
  layer = load(layer_fn);
  new_surface_gps_time = layer.GPS_time;
  new_surface = layer.layerData{1}.value{2}.data;
  % Get the previous frame if necessary
  if orig_records.gps_time(1) < new_surface_gps_time(1)-1
    layer_fn = fullfile(layer_path,sprintf('Data_%s_%03d.mat',param.day_seg,param.proc.frm-1));
    if exist(layer_fn,'file')
      layer = load(layer_fn);
      new_surface_gps_time = [layer.GPS_time new_surface_gps_time];
      new_surface = [layer.layerData{1}.value{2}.data new_surface];
    end
  end
  % Get the next frame if necessary
  if orig_records.gps_time(end) > new_surface_gps_time(end)+1
    layer_fn = fullfile(layer_path,sprintf('Data_%s_%03d.mat',param.day_seg,param.proc.frm+1));
    if exist(layer_fn,'file')
      layer = load(layer_fn);
      new_surface_gps_time = [new_surface_gps_time layer.GPS_time];
      new_surface = [new_surface layer.layerData{1}.value{2}.data];
    end
  end
  % Since layer files may have overlapping data, sort it
  [new_surface_gps_time new_surface_idxs] = sort(new_surface_gps_time);
  new_surface = new_surface(new_surface_idxs);
  
  % Do the interpolation and overwrite the orig_records.surface variable
  new_surface = interp1(new_surface_gps_time,new_surface,orig_records.gps_time,'linear','extrap');
  orig_records.surface = new_surface;
end

%% Load waveforms
% =========================================================================
if strcmpi(radar_name,'mcrds')
  [wfs,rec_data_size] = load_mcrds_wfs(orig_records.settings, param, ...
    1:max(old_param_records.records.file.adcs), param.csarp);
  load_param.load.rec_data_size = rec_data_size;
elseif any(strcmpi(radar_name,{'acords','hfrds','mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
  [wfs,rec_data_size] = load_mcords_wfs(orig_records.settings, param, ...
    1:max(old_param_records.records.file.adcs), param.csarp);
  load_param.load.rec_data_size = rec_data_size;
elseif any(strcmpi(radar_name,{'icards'}))% add icards----qishi
  [wfs,rec_data_size] = load_icards_wfs(orig_records.settings, param, ...
    1:max(old_param_records.records.file.adcs), param.csarp);
    load_param.load.rec_data_size = rec_data_size;
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5'}))
  wfs = load_fmcw_wfs(orig_records.settings, param, ...
    1:max(old_param_records.records.file.adcs), param.csarp);
end
load_param.wfs                = wfs;

%% Collect record file information required for using load_mcords_data
%  - Performs mapping between param.adcs and the records file contents
%  - Translates filenames from relative to absolute
%  - Makes filenames have the correct filesep
% =====================================================================

% Create a list of unique receivers required by the imgs list
param.load.adcs = [];
for idx = 1:length(param.load.imgs)
  param.load.adcs = unique(cat(2, ...
    abs(param.load.imgs{idx}(:,2)).', param.load.adcs));
end

recs = param.load.recs - param.load.recs(1) + 1;
if any(strcmpi(radar_name,{'hfrds','icards','mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
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
      param.load.recs,orig_records.relative_rec_num{board_idx});
    load_param.load.offset{idx} = orig_records.offset(board_idx,:);
    file_idxs = unique(load_param.load.file_idx{idx});
    
    % Recognize if first record is really from previous file and it is a
    % valid record (i.e. offset does not equal -2^31)
    if sign(load_param.load.offset{idx}(1)) < 0 && load_param.load.offset{idx}(1) ~= -2^31
      file_idxs = [file_idxs(1)-1 file_idxs];
    end
    
    % Just copy the filenames we need
    load_param.load.filenames{idx}(file_idxs) = orig_records.relative_filename{board_idx}(file_idxs);

    % Modify filename according to channel
    for file_idx = 1:length(load_param.load.filenames{idx})
      if any(strcmpi(radar_name,{'mcords5'}))
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
elseif strcmpi(radar_name,'mcrds')
  load_param.load.offset = orig_records.offset;
  load_param.load.file_rec_offset = orig_records.relative_rec_num;
  load_param.load.filenames = orig_records.relative_filename;
  base_dir = ct_filename_data(param,param.vectors.file.base_dir);
  adc_folder_name = param.vectors.file.adc_folder_name;
  load_param.load.filepath = fullfile(base_dir, adc_folder_name);
  load_param.load.wfs = orig_records.settings.wfs;
  load_param.load.wfs_records = orig_records.settings.wfs_records;  
elseif strcmpi(radar_name,'acords')
  load_param.load.offset = orig_records.offset;
  load_param.load.file_rec_offset = orig_records.relative_rec_num;
  load_param.load.filenames = orig_records.relative_filename;
  base_dir = ct_filename_data(param,param.vectors.file.base_dir);
  adc_folder_name = param.vectors.file.adc_folder_name;
  load_param.load.filepath = fullfile(base_dir, adc_folder_name);
  load_param.load.wfs = orig_records.settings.wfs;
  load_param.load.wfs_records = orig_records.settings.wfs_records;
elseif strcmpi(radar_name,'icards')% add icards---qishi
  load_param.load.offset = orig_records.offset;
  load_param.load.file_rec_offset = orig_records.relative_rec_num;
  load_param.load.filenames = orig_records.relative_filename;
  base_dir = ct_filename_data(param,param.vectors.file.base_dir);
  adc_folder_name = param.vectors.file.adc_folder_name;
  load_param.load.filepath = fullfile(base_dir, adc_folder_name);
  load_param.load.wfs = orig_records.settings.wfs;
  load_param.load.wfs_records = orig_records.settings.wfs_records; 
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5'}))
  % Determine which ADC boards are supported and which ones were actually loaded
  if ~isfield(old_param_records.records.file,'adc_headers') || isempty(old_param_records.records.file.adc_headers)
    old_param_records.records.file.adc_headers = old_param_records.records.file.adcs;
  end
  boards = adc_to_board(param.radar_name,old_param_records.records.file.adcs);
  boards_headers = adc_to_board(param.radar_name,old_param_records.records.file.adc_headers);
  
  for idx = 1:length(param.load.adcs)
    adc = param.load.adcs(idx);
    adc_idx = find(old_param_records.records.file.adcs == param.load.adcs(idx));
    if isempty(adc_idx)
      error('ADC %d not present in records file\n', param.load.adcs(idx));
    end
    
    % Just get the file-information for the records we need
    board = adc_to_board(param.radar_name,adc);
    actual_board_idx = find(board == boards);
    board_idx = find(old_param_records.records.file.adc_headers(actual_board_idx) == boards_headers);
    load_param.load.file_idx{idx} = relative_rec_num_to_file_idx_vector( ...
      load_param.load.recs,orig_records.relative_rec_num{board_idx});
    load_param.load.offset{idx} = orig_records.offset(board_idx,:);
    file_idxs = unique(load_param.load.file_idx{idx});
    
    % Recognize if first record is really from previous file and it is a
    % valid record (i.e. offset does not equal -2^31)
    if sign(load_param.load.offset{idx}(1)) < 0 && load_param.load.offset{idx}(1) ~= -2^31
      file_idxs = [file_idxs(1)-1 file_idxs];
    end
    
    % Just copy the filenames we need
    load_param.load.filenames{idx}(file_idxs) = orig_records.relative_filename{board_idx}(file_idxs);
    
    % Modify filename according to channel
    for file_idx = 1:length(load_param.load.filenames{idx})
      if ~isequal(old_param_records.records.file.adc_headers,old_param_records.records.file.adcs)
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

else
  error('Radar name %s not supported', param.radar_name);
end

%% Setup control parameters for loading data
% =====================================================================

load_param.load.adcs = param.load.adcs;

load_param.proc.pulse_comp         = param.csarp.pulse_comp;
load_param.proc.ft_dec             = param.csarp.ft_dec;
load_param.proc.ft_wind            = param.csarp.ft_wind;
load_param.proc.ft_wind_time       = param.csarp.ft_wind_time;
load_param.proc.presums            = param.csarp.presums;
load_param.proc.combine_rx         = param.csarp.combine_rx;
load_param.proc.pulse_rfi          = param.csarp.pulse_rfi;
load_param.proc.trim_vals          = param.csarp.trim_vals;
load_param.proc.coh_noise_method   = param.csarp.coh_noise_method;
load_param.proc.coh_noise_arg      = param.csarp.coh_noise_arg;

load_param.radar = param.radar;
load_param.surface = orig_records.surface;
if strcmpi(radar_name,'acords')
  load_param.load.file_version = param.records.file_version;
end

%% Load and Pulse Compress Data
% =====================================================================
% Load data into g_data using load_mcords_data
load_param.load.imgs = param.load.imgs;

if strcmpi(radar_name,'mcords')
  %   if strcmpi(param.season_name,'mcords_simulator')
  %     load_param.fn = get_filename(base_dir,'','','mat');
  %     load_simulated_data(load_param);
  %   else
  load_mcords_data(load_param);
  %   end
elseif any(strcmpi(radar_name,{'hfrds','mcords2','mcords3','mcords4','mcords5','seaice'}))
  load_mcords2_data(load_param);
elseif strcmpi(radar_name,'accum2')
  load_accum2_data(load_param);
elseif strcmpi(radar_name,'acords')
  load_acords_data(load_param);
elseif strcmpi(radar_name,'mcrds')
  if isfield(orig_records,'adc_phase_corr_deg') && isfield(param.radar,'adc_phase_corr_en') && param.radar.adc_phase_corr_en
    load_param.adc_phase_corr_deg = orig_records.adc_phase_corr_deg;
  else
    load_param.adc_phase_corr_deg = zeros(length(load_param.surface),max(orig_records.param_records.records.file.adcs));
  end
  load_mcrds_data(load_param);
elseif strcmpi(radar_name,'icards')
  load_icards_data(load_param,param);
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5'}))
  load_param.proc.elev_correction = 0;%param.csarp.elev_correction;
  load_param.proc.deconvolution = param.csarp.deconvolution;
  load_param.proc.psd_smooth = param.csarp.psd_smooth;
  load_param.radar_name = param.radar_name;
  load_param.season_name = param.season_name;
  load_param.day_seg = param.day_seg;
  load_param.load.tmp_path = param.tmp_path;
  load_param.out_path = param.out_path;
  [img_time,img_valid_rng,img_deconv_filter_idx,img_freq] = load_fmcw_data(load_param,orig_records);
  valid_rng = img_valid_rng{1};
  deconv_filter_idx = img_deconv_filter_idx{1};
  for wf = 1:length(wfs)
    wfs(wf).time = img_time{1};
    wfs(wf).freq = img_freq{1};
  end
end

%% Prepare trajectory information
% =========================================================================

%Decimate orig_records and ref according to presums
if param.csarp.presums > 1
  orig_records.lat = fir_dec(orig_records.lat,param.csarp.presums);
  orig_records.lon = fir_dec(orig_records.lon,param.csarp.presums);
  orig_records.elev = fir_dec(orig_records.elev,param.csarp.presums);
  orig_records.roll = fir_dec(orig_records.roll,param.csarp.presums);
  orig_records.pitch = fir_dec(orig_records.pitch,param.csarp.presums);
  orig_records.heading = fir_dec(orig_records.heading,param.csarp.presums);
  orig_records.gps_time = fir_dec(orig_records.gps_time,param.csarp.presums);
  orig_records.surface = fir_dec(orig_records.surface,param.csarp.presums);
  param.proc.along_track = fir_dec(param.proc.along_track,param.csarp.presums);
end

trajectory_param = struct('gps_source',orig_records.gps_source, ...
  'season_name',param.season_name,'radar_name',param.radar_name,'rx_path', 0, ...
  'tx_weights', [], 'lever_arm_fh', param.csarp.lever_arm_fh);
ref = trajectory_with_leverarm(orig_records,trajectory_param);

% Lsar = use approximate SAR aperture length
if isfield(param.csarp,'Lsar')
  Lsar = c/wfs(1).fc*(param.csarp.Lsar.agl+param.csarp.Lsar.thick/sqrt(er_ice))/(2*param.csarp.sigma_x);
else
  Lsar = c/wfs(1).fc*(500+1000/sqrt(er_ice))/(2*param.csarp.sigma_x);
end

along_track = param.proc.along_track - param.proc.along_track(1);

output_along_track = along_track(1) + param.proc.output_along_track_offset ...
  + param.csarp.sigma_x*(0:param.proc.output_along_track_Nx-1);

% Resample reference trajectory at output positions
% 1. Convert ref trajectory to ecef
ecef = zeros(3,size(ref.lat,2));
[ecef(1,:) ecef(2,:) ecef(3,:)] = geodetic2ecef(ref.lat/180*pi, ref.lon/180*pi, ref.elev, WGS84.ellipsoid);
% 2. Resample based on input and output along track
ecef = interp1(along_track,ecef.',output_along_track,'linear','extrap').';
% 3. Convert ecef to geodetic
[lat,lon,elev] = ecef2geodetic(ecef(1,:), ecef(2,:), ecef(3,:), WGS84.ellipsoid);
lat = lat*180/pi;
lon = lon*180/pi;

%% Remove coherent noise
% =========================================================================
if param.csarp.coh_noise_method && ~any(strcmpi(radar_name,{'kuband','snow','kuband2','snow2','kuband3','kaband3','snow3','snow5'}))
  
  if param.csarp.coh_noise_method == 3 && isempty(param.csarp.coh_noise_arg)
    param.csarp.coh_noise_arg = 255;
  end
  
  % DC and near-DC REMOVAL
  if param.csarp.ground_based
    % Only remove coherent noise from sections which are moving
    vel = diff(along_track) ./ diff(orig_records.gps_time);
    good_mask = vel > median(vel)/5;
    % figure(1); clf;
    % plot(vel);
    % hold on;
    % vel(~good_mask) = NaN;
    % plot(vel,'r');
    % hold off;
    for idx=1:length(g_data)
      % Transpose for faster memory access
      g_data{idx} = permute(g_data{idx},[2 1 3]);
      for rbin = 1:size(g_data{idx},2)
        for wf_adc_idx = 1:size(g_data{idx},3)
          if param.csarp.coh_noise_method == 1
            g_data{idx}(:,rbin,wf_adc_idx) = g_data{idx}(:,rbin,wf_adc_idx) - mean(g_data{idx}(good_mask,rbin,wf_adc_idx));
          else
            error('param.csarp.coh_noise_method %d not supported.',param.csarp.coh_noise_method);
          end
        end
      end
      % Undo transpose
      g_data{idx} = permute(g_data{idx},[2 1 3]);
    end
  else
    % Remove only the DC Doppler component
    for idx=1:length(g_data)
      for wf_adc_idx = 1:size(g_data{idx},3)
        if param.csarp.coh_noise_method == 1
          g_data{idx}(:,:,wf_adc_idx) = bsxfun(@minus, g_data{idx}(:,:,wf_adc_idx), ...
            mean(g_data{idx}(:,:,wf_adc_idx),2));
        elseif param.csarp.coh_noise_method == 3
          g_data{idx}(:,:,wf_adc_idx) = bsxfun(@minus, g_data{idx}(:,:,wf_adc_idx), ...
            fir_dec(g_data{idx}(:,:,wf_adc_idx),hanning(param.csarp.coh_noise_arg).'/(param.csarp.coh_noise_arg/2+0.5),1));
        else
          error('param.csarp.coh_noise_method %d not supported.',param.csarp.coh_noise_method);
        end
      end
    end
  end
end

%% Main loop to process each image
% =========================================================================
for img_idx = 1:length(load_param.load.imgs)
  if param.csarp.combine_rx
    % Receivers combined, so just get the first wf/adc pair info to name the outfile
    imgs_list = load_param.load.imgs{1}(1,:);
  else
    % Receivers processed individually, so get information for all wf/adc pairs.
    imgs_list = load_param.load.imgs{img_idx};
  end
  
  for wf_adc_idx = 1:size(imgs_list,1)
    % Processing loop
    % Runs once for combine_rx = true
    % Runs for each wf/adc pair in image if combine_rx = false
    
    wf = abs(imgs_list(wf_adc_idx,1));
    adc = abs(imgs_list(wf_adc_idx,2));
    
    %% Compute trajectory, SAR phase center and coordinate system
    if isempty(param.csarp.lever_arm_fh)
      records = orig_records;
    else
      if ~param.csarp.combine_rx
        % Create actual trajectory
        trajectory_param = struct('gps_source',orig_records.gps_source, ...
          'season_name',param.season_name,'radar_name',param.radar_name, ...
          'rx_path', wfs(wf).rx_paths(adc), ...
          'tx_weights', wfs(wf).tx_weights, 'lever_arm_fh', param.csarp.lever_arm_fh);
        records = trajectory_with_leverarm(orig_records,trajectory_param);
      else
        trajectory_param = struct('gps_source',orig_records.gps_source, ...
          'season_name',param.season_name,'radar_name',param.radar_name, ...
          'rx_path', wfs(wf).rx_paths(adc), ...
          'tx_weights', wfs(wf).tx_weights, 'lever_arm_fh', param.csarp.lever_arm_fh);
        for tmp_wf_adc_idx = 2:size(load_param.load.imgs{1},1)
          tmp_wf = abs(load_param.load.imgs{1}(tmp_wf_adc_idx,1));
          tmp_adc = abs(load_param.load.imgs{1}(tmp_wf_adc_idx,2));
          trajectory_param.rx_path(tmp_wf_adc_idx) = wfs(tmp_wf).rx_paths(tmp_adc);
        end
        records = trajectory_with_leverarm(orig_records,trajectory_param);
      end
    end
    
    SAR_coord_param.type = param.csarp.mocomp.type;
    SAR_coord_param.squint = [0 0 -1].';
    SAR_coord_param.Lsar = Lsar;
    fcs = SAR_coord_system(SAR_coord_param,records,ref,along_track,output_along_track);
    fcs.gps_source = orig_records.gps_source;
    
    %% SAR Processing
    if strcmpi(param.csarp.sar_type,'f-k')
      % f-k migration overview
      %
      % 1. Loop for each subaperture (repeat steps 2-4 for each aperture)
      %
      % 2. Motion compensation required before taking FFT
      %    g_data (raw with coherent noise optionally removed)
      %      --> g_data (motion compensated ft-fft)
      %
      % 3. Uniform re-sampling
      %   a. uniform_en = false, Assume data is uniformly sampled, apply fft in slow time and
      %      decimate in this domain by dropping doppler bins outside window
      %      g_data (motion compensated ft-fft)
      %        --> g_data (slow time ft/st-fft) [only done on the first loop]
      %        --> data (decimated ft/st-fft)
      %   b. uniform_en = true, Spatial filter to decimated axis and then take fft
      %      g_data (motion compensated ft-fft)
      %        --> data (decimated ft-fft)
      %        --> data (decimated ft/st-fft)
      %
      % 4. Regular f-k migration for each subaperture
      
      g_data{img_idx}(:,:,wf_adc_idx) = fft(g_data{img_idx}(:,:,wf_adc_idx),[],1);
      
      fcs.squint = [0 0 -1].';
      %fcs.squint = fcs.squint ./ sqrt(dot(fcs.squint,fcs.squint));
      
      %% Motion Compensation for f-k migration
      if param.csarp.mocomp.en
        % Determine the required motion compensation (drange and dx)
        %  Positive drange means the the range will be made longer, time delay
        %  will be made larger, and phase will be more negative
        %  Positive dx means that the data will be shifted forward (i.e. it
        %  currently lags behind)
        fcs.type = param.csarp.mocomp.type;
        fcs.filter = param.csarp.mocomp.filter;
        [drange,dx] = motion_comp(fcs,records,ref,along_track,output_along_track);
        
        % Time shift data in the frequency domain
        dtime = 2*drange/c;
        for rline = 1:size(g_data{img_idx},2)
          g_data{img_idx}(:,rline,wf_adc_idx) = g_data{img_idx}(:,rline,wf_adc_idx) ...
            .*exp(-1i*2*pi*wfs(wf).freq*dtime(rline));
        end
        
        along_track_mc = along_track + dx;
      else
        along_track_mc = along_track;
      end
      
      % output_along_track: these are the output values from csarp
      % output_along_track_pre/post: these are the buffers to keep the
      %   data from wrapping around in slow-time (i.e. linear convolution
      %   vs. circular convolution)
      % BUG!: At the beginning and end of a segment there is no data and
      % the buffer is not added... need to fix.
      output_along_track_pre = fliplr(output_along_track(1)-param.csarp.sigma_x:-param.csarp.sigma_x:along_track(1));
      output_along_track_post = output_along_track(end)+param.csarp.sigma_x:param.csarp.sigma_x:along_track(end);
      output_along_track_full = [output_along_track_pre output_along_track output_along_track_post];

      %% Prepare subaperture variables
      num_subapertures = length(param.csarp.sub_aperture_steering);
      if mod(num_subapertures,2) ~= 1
        error('Number of subapertures must be even');
      end
      if any(param.csarp.sub_aperture_steering ~= -(num_subapertures-1)/4 : 0.5 : (num_subapertures-1)/4)
        error('param.csarp.sub_aperture_steering must be of the form -N:0.5:N');
      end
      proc_oversample = (1+num_subapertures)/2; % Oversampling rate
      proc_sigma_x = param.csarp.sigma_x / proc_oversample;
      proc_along_track = output_along_track_full(1) ...
        + proc_sigma_x * (0:length(output_along_track_full)*proc_oversample-1);

      %% Uniform resampling and subaperture selection for f-k migration
      if param.csarp.mocomp.uniform_en
        % Uniformly resample data in slow-time onto output along-track
        data = arbitrary_resample(g_data{img_idx}(:,:,wf_adc_idx), ...
          along_track_mc,proc_along_track, struct('filt_len', ...
          proc_sigma_x*16,'dx',proc_sigma_x,'method','sinc'));
        data = fft(data,[],2);
        
      else
        % Assume data is already uniformly sampled in slow-time
        % There are lots of approximations in this section... it's a bit
        % of a hack.
        x_lin = linspace(along_track(1), ...
          along_track(end),length(along_track));
        
        % Create kx (along-track wavenumber) axis of input data
        kx = gen_kx(x_lin);
        
        % Create kx axis of output (desired) data
        kx_desired = gen_kx(output_along_track_full);
        
        % Take slow-time FFT and decimate the data onto the desired
        % output along track positions by selecting just the doppler
        % bins that correspond to this

        g_data{img_idx}(:,:,wf_adc_idx) = fft(g_data{img_idx}(:,:,wf_adc_idx),[],2);
        filt_idx = kx < max(kx_desired) & kx > min(kx_desired);
        % Since kx and kx_desired won't match up perfectly, we may have
        % to append a few more doppler bins to get the numbers to line
        % up.
        if length(output_along_track_full) - sum(filt_idx) == 1
          filt_idx(find(filt_idx==0,1)) = 1;
        elseif length(output_along_track_full) - sum(filt_idx) == 2
          filt_idx(find(filt_idx==0,1)) = 1;
          filt_idx(find(filt_idx==0,1,'last')) = 1;
        end
        filt_len = length(find(filt_idx));
        filt_idx = filt_idx.';
        filt_idx = find(filt_idx);
        kx_sa    = kx(filt_idx);
        [kx_sa,kx_idxs] = sort(kx_sa);
        filt_idx = filt_idx(kx_idxs);
        filt_idx = ifftshift(filt_idx);
        data = g_data{img_idx}(:,filt_idx,wf_adc_idx);
      end
      
      %% F-k migration
      if param.csarp.mocomp.en
        eps_r  = perm_profile(mean(records.surface + dtime),wfs(wf).time,'constant',param.csarp.start_eps);
      else
        eps_r  = perm_profile(mean(records.surface),wfs(wf).time,'constant',param.csarp.start_eps);
      end

      kx = gen_kx(proc_along_track);
      fk_data_ml = fk_migration(data,wfs(wf).time,wfs(wf).freq,kx,eps_r,param);
      fk_data_ml = fk_data_ml(:,1+length(output_along_track_pre):end-length(output_along_track_post),:);

      if 0
        % DEBUG code
        figure(1); clf;
        for subap = 1:num_subapertures
          imagesc(lp(fk_data_ml(300:700,:,subap)))
          title(sprintf('%d',subap ));
          pause(1);
        end
        figure(1); clf;
        imagesc(lp(mean(abs(fk_data_ml(300:700,:,:)).^2,3)));
        figure(2); clf;
        imagesc(lp(mean(abs(fk_data_ml(300:700,:,1:6)).^2,3)));
        figure(3); clf;
        imagesc(lp(mean(abs(fk_data_ml(300:700,:,end-5:end)).^2,3)));
      end
      
      if param.csarp.mocomp.en
        %% Undo motion compensation
        % Resample dtime to f-k migration output
        dtime = interp1(orig_records.gps_time,dtime,fcs.gps_time);
        
        % Time shift data in the frequency domain
        fk_data_ml = fft(fk_data_ml,[],1);
        fk_data_ml = fk_data_ml.*exp(1i*2*pi*repmat(wfs(wf).freq, [1,size(fk_data_ml,2),size(fk_data_ml,3)]) ...
          .*repmat(dtime, [size(fk_data_ml,1),1,size(fk_data_ml,3)]));
        fk_data_ml = ifft(fk_data_ml,[],1);
      end
      
      %% Save Radar data
      for subap = 1:num_subapertures
        % Create output path
        out_path = fullfile(ct_filename_out(param, ...
          param.csarp.out_path, 'CSARP_out'), ...
          sprintf('fk_data_%03d_%02d_%02d',param.proc.frm, ...
          subap, param.proc.sub_band_idx));
        if ~exist(out_path,'dir')
          mkdir(out_path);
        end
        
        % Create filename
        % - Hack: multiple receivers are named with the first receiver in the list
        out_fn = sprintf('wf_%02d_adc_%02d_chk_%03d',wf, adc,param.csarp.chunk_id);
        out_full_fn = fullfile(out_path,[out_fn '.mat']);
        
        % Save
        fprintf('  Saving output %s\n', out_full_fn);
        param_records = old_param_records;
        param_csarp = param;
        fk_data = fk_data_ml(:,:,subap);
        save('-v7.3',out_full_fn,'fk_data','fcs','lat','lon','elev','wfs','param_csarp','param_records');
      end
      
    elseif strcmpi(param.csarp.sar_type,'tdbp')
      % time domain backporjection overview
      data = g_data{img_idx}(:,:,wf_adc_idx);
      
      % set up SAR coordinate system
      [x_ecef, y_ecef, z_ecef] = geodetic2ecef(records.lat*pi/180, records.lon*pi/180, records.elev, WGS84.ellipsoid);
      records.lon_ref = mean(records.lon);
      records.lat_ref = mean(records.lat);
      records.elev_ref = mean(records.elev);
      [x_enu,y_enu,z_enu] = ecef2lv(x_ecef, y_ecef, z_ecef, records.lat_ref*pi/180, records.lon_ref*pi/180, records.elev_ref, WGS84.ellipsoid);
%       along_track =  [0 cumsum(sqrt(diff(x_enu).^2 + diff(y_enu).^2))];
      SAR_coord_param.phase_center = [x_enu;y_enu;z_enu];
      SAR_coord_param.Lsar = Lsar;
      SAR_coord_param.along_track = along_track;
      SAR_coord_param.output_along_track_idxs = param.proc.output_along_track_idxs;
%       if strcmpi(param.season_name,'mcords_simulator') % for 20110708_01_001 simulated data
%         param.proc.output_along_track_idxs = [fliplr([4632:-10:0]),[4642:10:length(along_track)]];
%         SAR_coord_param.output_along_track_idxs = param.proc.output_along_track_idxs;
%       end
      SAR_coord_param.wfs = wfs(wf);
      
      % surface tracker
      % two methods to get ice surface: param.surf_source 1/2
      % 1: from get_heights; 2:from laser data;
      param.surf_source = 1;
      if param.surf_source == 1
        surfTimes = records.surface;
      elseif param.surf_source == 2
        param.laser_surface = 1;
        param.laser_data_fn = '/cresis/projects/metadata/2008_Greenland_TO_icessn/2008_Greenland/080801a_icessn_nadir0seg';
        param.laser_data_fn = '/cresis/projects/metadata/2008_Greenland_TO_icessn/2008_Greenland/080707_icessn_nadir0seg';
        fid = fopen(param.laser_data_fn);
        [laser_data_tmp] = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
        fclose(fid);
        Year = 2008;
        Mon = 7;
        Day = 7;
        laser_data.gps_time = (datenum(Year,Mon,Day)-datenum(1970,1,1))*86400 + laser_data_tmp{1};
        laser_data.surf_elev = laser_data_tmp{4};
        laser_data.surf_elev = interp1(laser_data.gps_time,laser_data.surf_elev,records.gps_time);
        surfTimes = 2*(records.elev-laser_data.surf_elev)/c;
        clear laser_data_tmp;
      end
      
      SAR_coord_param.surf = zeros(3,length(along_track));
      SAR_coord_param.surf(1,:) = x_enu;
      SAR_coord_param.surf(2,:) = y_enu;
      SAR_coord_param.surf(3,:) = z_enu - surfTimes*c/2;
      SAR_coord_param.surf_p = polyfit(along_track,SAR_coord_param.surf(3,:),10);
      SAR_coord_param.surf(3,:) = polyval(SAR_coord_param.surf_p,along_track);
      N = length(SAR_coord_param.surf_p);
      surfSlope = SAR_coord_param.surf_p(N-1);
      x_pwr = along_track;
      for ii = N-3:-1:1
        surfSlope = surfSlope + (N-ii)*SAR_coord_param.surf_p(ii)*along_track.*x_pwr;
        x_pwr = x_pwr.*along_track;
      end
      surfSlope = surfSlope + 2*SAR_coord_param.surf_p(N-2)*along_track;
      surfSlope(abs(surfSlope)<0.0001*pi/180) = 0;
      SAR_coord_param.surfAngle = atan(surfSlope);
      % SAR_coord_param.surfNormal = zeros(3,length(along_track));
      % SAR_coord_param.surfNormal(1,:) = cos(atan(surfSlope)+pi/2);
      % SAR_coord_param.surfNormal(3,:) = sin(atan(surfSlope)+pi/2);
      % SAR_coord_param.surfNormalAngle = atan(surfSlope)+pi/2;
      SAR_coord_param.surfBins = round((2*(z_enu-SAR_coord_param.surf(3,:))/c-wfs(wf).time(1))/wfs(wf).dt) + 1;
      
      n = size(data,1);
      m = length(param.proc.output_along_track_idxs);
      SAR_coord_param.pixel = zeros(3,n,m);
      SAR_coord_param.pixel(1,:,:) = repmat(x_enu(param.proc.output_along_track_idxs),n,1);
      SAR_coord_param.pixel(2,:,:) = repmat(y_enu(param.proc.output_along_track_idxs),n,1);
      eta_ice = sqrt(er_ice);
      for line = 1:m
        out_idx = param.proc.output_along_track_idxs(line);
        SAR_coord_param.pixel(3,1:SAR_coord_param.surfBins(out_idx),line) = z_enu(out_idx) - wfs(wf).time(1)*c/2 - ...
          [(0:SAR_coord_param.surfBins(out_idx)-1)'*wfs(wf).dt*c/2];
        SAR_coord_param.pixel(3,SAR_coord_param.surfBins(out_idx)+1,line) = ...
          SAR_coord_param.pixel(3,SAR_coord_param.surfBins(out_idx),line) -...
          c*(surfTimes(out_idx)-wfs(wf).time(SAR_coord_param.surfBins(out_idx)))/2 -...
          (wfs(wf).time(SAR_coord_param.surfBins(out_idx)+1)-surfTimes(out_idx))*c/(2*eta_ice);
        SAR_coord_param.pixel(3,SAR_coord_param.surfBins(out_idx)+2:n,line) = ...
          SAR_coord_param.pixel(3,SAR_coord_param.surfBins(out_idx)+1,line) - ...
          (1:n-SAR_coord_param.surfBins(out_idx)-1)'*wfs(wf).dt*c/(2*eta_ice);
      end
      SAR_coord_param.h = SAR_coord_param.phase_center(3,:)-SAR_coord_param.surf(3,:);
      SAR_coord_param.h_mean = mean(SAR_coord_param.h);
      Lsar_surf = c/wfs(wf).fc*SAR_coord_param.h_mean/(2*param.csarp.sigma_x);
      SAR_coord_param.HbeamWidth = atan(0.5*Lsar_surf/SAR_coord_param.h_mean);
      
      tdbp_param = SAR_coord_param;
      clear SAR_coord_param;
      tdbp_param.proc.Nfft = 2^ceil(log2(length(wfs(wf).time_raw)));
      tdbp_param.proc.skip_surf = param.csarp.skip_surf;
      tdbp_param.proc.start_range_bin_above_surf = param.csarp.start_range_bin_above_surf;
      tdbp_param.proc.start_range_bin = param.csarp.start_range_bin;
      tdbp_param.proc.end_range_bin = param.csarp.end_range_bin;
      tdbp_param.refraction_method = param.csarp.refraction_method;
      tdbp_param.fc = wfs(wf).fc;
      tdbp_param.time = wfs(wf).time;
      tdbp_param.c = c;
      tdbp_param.eta_ice = eta_ice;
      tdbp_param.st_wind = param.csarp.st_wind;
      tdbp_param.sub_aperture_steering = param.csarp.sub_aperture_steering;
            
      fcs.squint = [0 0 -1].';
      tdbp_data0 = tdbp(tdbp_param,data);
      
      for subap = 1:size(tdbp_data0,3) % save each subaperture data to its own folder
        % Create output path
        out_path = fullfile(ct_filename_out(param,param.csarp.out_path, 'CSARP_out'),...
          sprintf('tdbp_data_%03d_%02d_%02d',param.proc.frm,subap, param.proc.sub_band_idx));
        if ~exist(out_path,'dir')
          mkdir(out_path);
        end
        
        % Create filename
        % - Hack: multiple receivers are named with the first receiver in the list
        out_fn = sprintf('wf_%02d_adc_%02d_chk_%03d', wf, adc, param.csarp.chunk_id);
        out_full_fn = fullfile(out_path,[out_fn '.mat']);
        
        % Save
        fprintf('  Saving output %s\n', out_full_fn);
        param_records = old_param_records;
        param_csarp = param;
        param_csarp.tdbp = tdbp_param;
        tdbp_data = tdbp_data0(:,:,subap);
        save('-v7.3',out_full_fn,'tdbp_data','fcs','lat','lon','elev','wfs','param_csarp','param_records','tdbp_param');
      end
    elseif strcmpi(param.csarp.sar_type,'mltdp')
      fcs.squint = [0 0 -1].';
      [B,A] = butter(4,0.1);

      % Force elevation to be smooth (might be required for refraction)
      smoothed_elevation = filtfilt(B,A,records.elev);
      smoothed_ref_elevation = filtfilt(B,A,ref.elev);

      % Fit surface to polynomial to force it to be smooth (required for refraction)
      %  - Fit is done with special x-axis to prevent bad conditioning
      smoothed_surface = filtfilt(B,A,records.surface);
      sz_poly_order = 11;
      xfit = linspace(-1,1,length(smoothed_surface));
      smoothed_surface = polyval(polyfit(xfit,smoothed_surface,sz_poly_order),xfit);
      if 0 % set to 1 for surface fit over whole frame
        smoothed_elevation = interp1(param.proc.along_track_frm,param.proc.smoothed_elevation,param.proc.along_track,'linear','extrap');
        smoothed_surface = interp1(param.proc.along_track_frm,param.proc.smoothed_surface,param.proc.along_track,'linear','extrap');
      end

      data = g_data{img_idx}(:,:,wf_adc_idx);        
      % options for processing window in fast time
      surfBins_at_output = round((smoothed_surface-wfs(wf).time(1))/wfs(wf).dt)+1;
      if isempty(param.csarp.skip_surf)
        param.scarp.skip_surf = 0;                   % default value
      end
      if isempty(param.csarp.start_range_bin_above_surf)
        param.csarp.start_range_bin_above_surf = 5;  % default value
      end
      if param.csarp.skip_surf
        if isempty(param.csarp.start_range_bin)
          param.csarp.start_range_bin = max(surfBins_at_output) + 5;   % default value
        end
      else
        param.csarp.start_range_bin = min(surfBins_at_output) - param.csarp.start_range_bin_above_surf;        % default value
      end
      if isempty(param.csarp.end_range_bin)
        param.csarp.end_range_bin = size(data,1);  % default value
      end
%       if strcmpi(param.season_name,'mcords_simulator') % for 20110708_01_001 simulated data
%         output_along_track(463) = along_track(4632);
%         output_along_track(1:462) = output_along_track(463)-[462:-1:1]*param.csarp.sigma_x;
%         output_along_track(464:925) = output_along_track(463)+[1:462]*param.csarp.sigma_x;
%       end
      mltdp_data0 = ml_tdp(data,wfs(wf).fc,wfs(wf).time,along_track, smoothed_ref_elevation, ...
        smoothed_elevation,smoothed_surface,output_along_track,Lsar, ...
        length(param.csarp.sub_aperture_steering),param.csarp.start_range_bin,param.csarp.end_range_bin,param.csarp.start_eps);
      
      for subap = 1:size(mltdp_data0,3) % save each subaperture data to its own folder
        % Create output path
        out_path = fullfile(ct_filename_out(param,param.csarp.out_path, 'CSARP_out'),...
          sprintf('mltdp_data_%03d_%02d_%02d',param.proc.frm,subap, param.proc.sub_band_idx));
        if ~exist(out_path,'dir')
          mkdir(out_path);
        end
        
        % Create filename
        % - Hack: multiple receivers are named with the first receiver in the list
        out_fn = sprintf('wf_%02d_adc_%02d_chk_%03d',wf, adc,param.csarp.chunk_id);
        out_full_fn = fullfile(out_path,[out_fn '.mat']);
        
        fprintf('  Saving output %s\n', out_full_fn);
        param_records = old_param_records;
        param_csarp = param;
        mltdp_data = mltdp_data0(:,:,subap);
        save('-v7.3',out_full_fn,'mltdp_data','fcs','lat','lon','elev','wfs','param_csarp','param_records');
      end
    end
  end
end
clear data_storage;
success = true;
return
