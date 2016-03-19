function [success] = coh_noise_tracker_task(param)
% [success] = coh_noise_tracker_task(param)
%
% Cluster task for coh_noise_tracker. Does the actual data loading
% and noise analysis.
%
% param = struct controlling the loading, processing, surface tracking,
%   and quick look generation
%  .load = structure for which records to load
%   .records_fn = filename of records file
%   .recs = current records
%   .imgs = cell vector of images to load, each image is Nx2 array of
%     wf/adc pairs
%     NOTE: wfs/adc pairs are not indices into anything, they are the absolute
%     waveform/adc numbers. The records file will be loaded to decode
%     which index each wf/adc belongs to.
%  .debug_level = debug level (scalar integer)
%
%  .proc = structure containing information about framing
%   .frm = only used to determine the filename of the output
%
% .get_heights = structure controlling get_heights processing
%  .radar_name = name of radar string
%  .season_name = name of mission string
%  .day_seg = day-segment string
%  
%  get_heights fields used by load_mcords_data.m (see that function for details)
%  .ft_wind
%  .ft_wind_time
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
%  get_heights fields for post processing
%  .roll_correction = boolean, whether or not to apply roll phase correction
%  .lever_arm_fh = string containing function name
%  .elev_correction = boolean, whether or not to apply elevation phase correction
%  .B_filter = double vector, FIR filter coefficients to apply before
%    decimating, this function loads data before and after this frame
%    (if available) to avoid transients at the beginning and end
%  .decimate_factor = positive integer, decimation rate
%  .inc_ave = positive integer, number of incoherent averages to apply
%    (also decimates by this number)
%
%  .surf = get_heights structure controlling surface tracking
%   .en = boolean, whether or not to apply surface tracking
%   .wf_idx = positive integer, which waveform in the wfs list to apply surface tracking on
%   .min_bin = double scalar, the minimum range time that the surface can be tracked to.
%     This is used to keep the surface tracking routine from picking up the
%     feedthrough.  It requires a minimum elevation AGL.
%   .manual = boolean, whether or not to enable the manual tracking
%     interface.  Generally better to let the automated routine run, fix in
%     picker, and then update records (so surf.manual is mostly for debugging)
%
%  .qlook = get_heights structure controlling quick look generation
%   .en = boolean, whether or not to produce a quick look product
%    .out_path = output path of the quick look.  Three forms:
%      1. empty: default path based on the gRadar.out_path, param.records_name,
%         param.season_name, and param.day_seg
%      2. relative path: path based on gRadar.out_path and the contents of
%         .qlook.out_path
%      3. absolute path: uses this specific path for outputs
%   .wf_comb = vector of times of when to combine waveforms
%
% success = boolean which is true when the function executes properly
%   if a task fails before it can return success, then success will be
%   empty
% surfTimes = vector of propagation delays to the surface
%
% Author: John Paden
%
% See also coh_noise_tracker.m

global g_data;

physical_constants;
surfTimes = [];

records_fn = ct_filename_support(param,param.records.records_fn,'records');

if ~isfield(param.get_heights,'elev_correction') || isempty(param.get_heights.elev_correction)
  param.get_heights.elev_correction = false;
end

if ~param.get_heights.elev_correction ...
    && length(param.get_heights.B_filter) == param.get_heights.decimate_factor ...
    && all(param.get_heights.B_filter == 1)
  simple_firdec = true;
else
  simple_firdec = false;
end

if ~isfield(param.get_heights,'trim_vals') || isempty(param.get_heights.trim_vals)
  param.get_heights.trim_vals = [0 0];
end

if ~isfield(param.get_heights,'coh_noise_method') || isempty(param.get_heights.coh_noise_method)
  param.get_heights.coh_noise_method = 0;
end

if ~isfield(param.get_heights,'pulse_rfi')
  param.get_heights.pulse_rfi.en = 0;
end

if ~isfield(param.get_heights,'ft_dec')
  param.get_heights.ft_dec = 1;
end

if ~isfield(param.records,'file_version')
  param.records.file_version = [];
end

% =====================================================================
% Determine which records to load with load_mcords_data
%
% Load records on either side of the current block, note if at the
% beginning or end of the segment.  Load with minimal presumming.

load_param.load.recs(1) = param.load.recs(1);
load_param.load.recs(2) = param.load.recs(2);
records = read_records_aux_files(records_fn,load_param.load.recs);
old_param_records = records.param_records;
param_records = records.param_records;
param_records.gps_source = records.gps_source;

% =====================================================================
% Collect waveform information into one structure
%  (used by load_RADARNAME_data)

if simple_firdec
  param.get_heights.presums = param.get_heights.decimate_factor;
else
  param.get_heights.presums = 1;
end
param.get_heights.pulse_comp = 1;
param.get_heights.combine_rx = 0;
if strcmpi(param.radar_name,'mcrds')
  [wfs,rec_data_size] = load_mcrds_wfs(records.settings, param, ...
    1:max(param_records.records.file.adcs), param.get_heights);
  load_param.load.rec_data_size = rec_data_size;
elseif any(strcmpi(param.radar_name,{'mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
  [wfs,rec_data_size] = load_mcords_wfs(records.settings, param, ...
    1:max(param_records.records.file.adcs), param.get_heights);
  load_param.load.rec_data_size = rec_data_size;
elseif any(strcmpi(param.radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3'}))
  [path name] = fileparts(records_fn);
  cdf_fn = fullfile(path, sprintf('%s.nc', name));
  try
    records.settings.nyquist_zone = ncread(cdf_fn,'settings(1).nyquist_zone', ...
      [1 load_param.load.recs(1)],[1 diff(load_param.load.recs([1 end]))+1]);
  end
  try
    records.settings.loopback_mode = ncread(cdf_fn,'settings(1).loopback_mode',[1 load_param.load.recs(1)],[1 1]);
  end
  wfs_idx = find(records.settings.wfs_records <= load_param.load.recs(1),1,'last');
  records.settings.wfs = records.settings.wfs(wfs_idx).wfs;
  wfs = load_fmcw_wfs(records.settings, param, ...
    1:max(param_records.records.file.adcs), param.get_heights);
end
% load_wf_cmd=['load_',param.radar_name,'_wfs(records.wfs,param,1:max(param_records.file.adcs),param.get_heights)'];
% [wfs,rec_data_size] = eval(load_wf_cmd);
load_param.wfs = wfs;

% =====================================================================
% Collect record file information required for using load_RADAR_NAME_data
%  - Performs mapping between param.rxs and the records file contents
%  - Translates filenames from relative to absolute
%  - Makes filenames have the correct filesep
 
% Create a list of unique adcs required by the imgs list
param.load.adcs = [];
for idx = 1:length(param.load.imgs)
  new_adcs = abs(param.load.imgs{idx}(:,2:2:end)).';
  param.load.adcs = unique(cat(2, new_adcs(:).', param.load.adcs));
end

recs = load_param.load.recs - load_param.load.recs(1) + 1;
if any(strcmpi(param.radar_name,{'mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
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
    load_param.load.offset{idx} = records.offset(board_idx,:);
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
  load_param.load.offset = records.offset;
  load_param.load.file_rec_offset = records.relative_rec_num;
  load_param.load.filenames = records.relative_filename;
  base_dir = ct_filename_data(param,param.vectors.file.base_dir);
  adc_folder_name = param.vectors.file.adc_folder_name;
  load_param.load.filepath = fullfile(base_dir, adc_folder_name);
  load_param.load.wfs = records.settings.wfs;
  load_param.load.wfs_records = records.settings.wfs_records;
elseif any(strcmpi(param.radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3'}))
  load_param.load.offset{1} = records.offset;
  load_param.load.file_rec_offset = records.relative_rec_num;
  load_param.load.filenames = records.relative_filename;
  base_dir = ct_filename_data(param,param.vectors.file.base_dir);
  adc_folder_name = param.vectors.file.adc_folder_name;
  load_param.load.filepath = fullfile(base_dir, adc_folder_name);
  load_param.load.wfs = records.settings.wfs;
  load_param.load.radar_name = param.radar_name;
  load_param.load.season_name = param.season_name;
  load_param.load.day_seg = param.day_seg;
  load_param.load.tmp_path = param.tmp_path;
  load_param.load.file_version = param.records.file_version;
  load_param.load.file_idx{idx} = relative_rec_num_to_file_idx_vector( ...
    load_param.load.recs,records.relative_rec_num{1});
  
  filepath = get_segment_file_list(param,1);
  
  % Convert relative file paths into absolute file paths if required,
  % also corrects filesep (\ and /)
  for file_idx = 1:length(load_param.load.filenames{idx})
    load_param.load.filenames{idx}{file_idx} ...
      = fullfile(filepath,load_param.load.filenames{idx}{file_idx});
  end
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
else
  error('Radar name %s not supported', param.radar_name);
end

% =====================================================================
% Setup control parameters for load_mcords_data

load_param.load.adcs = param.load.adcs;

load_param.proc.trim_vals           = param.get_heights.trim_vals;
load_param.proc.pulse_comp          = param.get_heights.pulse_comp;
load_param.proc.ft_dec              = param.get_heights.ft_dec;
load_param.proc.ft_wind             = param.get_heights.ft_wind;
load_param.proc.ft_wind_time        = param.get_heights.ft_wind_time;
load_param.proc.presums             = param.get_heights.presums;
load_param.proc.combine_rx          = param.get_heights.combine_rx;
load_param.proc.pulse_rfi           = param.get_heights.pulse_rfi;
load_param.proc.coh_noise_method    = param.get_heights.coh_noise_method;
load_param.proc.coh_noise_arg       = param.get_heights.coh_noise_arg;

if param.analysis.coh_ave.en || param.analysis.specular.en
  load_param.proc.presums             = 1;
  param.get_heights.decimate_factor   = 1;
end

%   load_param.proc.combine_rx          = param.get_heights.combine_rx;
load_param.proc.pulse_rfi.en        = 0;

load_param.radar = param.radar;
load_param.surface = records.surface;

% =====================================================================
% Load and process each waveform separately
%
% For each waveform:
% 1. Load receiver data separately (minimal presumming)
% 2. Remove coherent noise (slow-time mean removal)
% 3. Apply roll correction
% 4. Combine receivers
% 5. Apply elevation compensation
% 6. FIR decimate the data
% =====================================================================

for img_idx = 1:length(param.load.imgs)
  
  % Default values to use
  wf = abs(param.load.imgs{img_idx}(1,1));
  adc = abs(param.load.imgs{img_idx}(1,2));
  lambda_fc = c/wfs(wf).fc;
  
  % Apply lever arm correction to trajectory data
  trajectory_param = struct('gps_source',param_records.gps_source, ...
    'season_name',param.season_name,'radar_name',param.radar_name, ...
    'rx_path', wfs(wf).rx_paths(adc), ...
    'tx_weights', wfs(wf).tx_weights, 'lever_arm_fh', param.get_heights.lever_arm_fh);
  for tmp_wf_adc_idx = 2:size(param.load.imgs{1},1)
    tmp_wf = abs(param.load.imgs{1}(tmp_wf_adc_idx,1));
    tmp_adc = abs(param.load.imgs{1}(tmp_wf_adc_idx,2));
    trajectory_param.rx_path(tmp_wf_adc_idx) = wfs(tmp_wf).rx_paths(tmp_adc);
  end
  out_records = trajectory_with_leverarm(records,trajectory_param);

  %% Load data into g_data using load_mcords_data
  load_param.load.imgs = param.load.imgs(img_idx);
  % Determine combination times when multiple wf-adc pairs are being loaded
  % to form a single range line
  if size(load_param.load.imgs{1},2) == 2
    load_param.load.wf_adc_comb.en = 0;
  else
    load_param.load.wf_adc_comb.en = 1;
    wf1 = load_param.load.imgs{1}(1,1);
    wf2 = load_param.load.imgs{1}(1,3);
    % t1 = time to switch from wf1 to wf2
    %wf_adc_surface = fir_dec(records.surface, param.get_heights.decimate_factor);
    wf_adc_surface = records.surface;
    t1 = eval(param.get_heights.qlook.wf_adc_comb{1});
    % If the time to switch is longer than wf1 then it gets capped to wf1
    t1(t1 > wfs(wf1).time(end) - wfs(wf1).Tpd) = wfs(wf1).time(end) - wfs(wf1).Tpd;
    load_param.load.wf_adc_comb.Nt = round((wfs(wf2).time(end) - wfs(wf1).time(1)) / wfs(wf1).dt);
    load_param.load.wf_adc_comb.rbins(1,:) = round((t1-wfs(wf1).time(1))/wfs(wf1).dt);
    load_param.load.wf_adc_comb.rbins(2,:) = 1+round((t1-wfs(wf2).time(1))/wfs(wf2).dt);
    % 1:load_param.load.img_comb{1}(1,rline) <- 1:load_param.load.img_comb{1}(1,rline)
    % load_param.load.img_comb{1}(1,rline)+1:end <- load_param.load.img_comb{1}(2,rline):end
  end
  if strcmpi(param.radar_name,'mcords')
    load_mcords_data(load_param);
    g_data = g_data{1};
  elseif any(strcmpi(param.radar_name,{'mcords2','mcords3','mcords4','mcords5'}))
    load_mcords2_data(load_param);
    g_data = g_data{1};
  elseif strcmpi(param.radar_name,'mcrds')
    if isfield(records,'adc_phase_corr_deg') && isfield(param.radar,'adc_phase_corr_en') && param.radar.adc_phase_corr_en
      load_param.adc_phase_corr_deg = records.adc_phase_corr_deg;
    else
      load_param.adc_phase_corr_deg = zeros(length(load_param.surface),max(records.param_records.records.file.adcs));
    end
    load_mcrds_data(load_param);
    g_data = g_data{1};
  elseif strcmpi(param.radar_name,'accum2')
    load_accum2_data(load_param);
    g_data = g_data{1};
  elseif any(strcmpi(param.radar_name,{'kuband','snow','kuband2','snow2','kuband3','snow3'}))
    load_param.proc.deconvolution = 0;
    load_param.proc.elev_correction = 0;
    load_param.proc.psd_smooth = 0;
    if param.analysis.coh_ave.en
      load_param.proc.pulse_comp = 0;
    end
    load_param.proc.coh_noise_tracker = 1;
    load_param.radar_name = param.radar_name;
    load_param.season_name = param.season_name;
    load_param.day_seg = param.day_seg;
    load_param.out_path = param.out_path;
    [img_time,img_valid_rng,img_deconv_filter_idx,img_freq,img_Mt] = load_fmcw_data(load_param,out_records);
    wfs.time = img_time;
  end
  
  if param.analysis.surf.en
    %% 1. Load layer
    layer_params = [];
    idx = 1;
    layer_params(idx).name = 'surface';
    layer_params(idx).source = 'records';
    layers = opsLoadLayers(param,layer_params);
    
    %% Apply lever arm correction to trajectory data, but preserve each
    %% channel separately.
    trajectory_param = struct('gps_source',param_records.gps_source, ...
      'season_name',param.season_name,'radar_name',param.radar_name, ...
      'rx_path', wfs(wf).rx_paths(adc), ...
      'tx_weights', wfs(wf).tx_weights, 'lever_arm_fh', param.get_heights.lever_arm_fh);
    out_records = trajectory_with_leverarm(records,trajectory_param);
    for tmp_wf_adc_idx = 2:size(param.load.imgs{1},1)
      tmp_wf = abs(param.load.imgs{1}(tmp_wf_adc_idx,1));
      tmp_adc = abs(param.load.imgs{1}(tmp_wf_adc_idx,2));
      trajectory_param.rx_path = wfs(tmp_wf).rx_paths(tmp_adc);
      tmp_records = trajectory_with_leverarm(records,trajectory_param);
      % Add the positions to the existing out_records
      out_records.gps_time = cat(1,out_records.gps_time,tmp_records.gps_time);
      out_records.lat = cat(1,out_records.lat,tmp_records.lat);
      out_records.lon = cat(1,out_records.lon,tmp_records.lon);
      out_records.elev = cat(1,out_records.elev,tmp_records.elev);
      out_records.roll = cat(1,out_records.roll,tmp_records.roll);
      out_records.pitch = cat(1,out_records.pitch,tmp_records.pitch);
      out_records.heading = cat(1,out_records.heading,tmp_records.heading);
    end
    
    %% FIR Decimate
    if simple_firdec
      out_records.gps_time = fir_dec(out_records.gps_time, param.get_heights.decimate_factor);
      out_records.lat = fir_dec(out_records.lat, param.get_heights.decimate_factor);
      out_records.lon = fir_dec(out_records.lon, param.get_heights.decimate_factor);
      out_records.elev = fir_dec(out_records.elev, param.get_heights.decimate_factor);
      out_records.roll = fir_dec(out_records.roll, param.get_heights.decimate_factor);
      out_records.pitch = fir_dec(out_records.pitch, param.get_heights.decimate_factor);
      out_records.heading = fir_dec(out_records.heading, param.get_heights.decimate_factor);
      
    else
      % Check for edge conditions that caused not enough data to be loaded
      % in the case where a segment starts and stops.
      % If not enough data was loaded, modify the filter coefficient
      % normalization so that it handles this
      Nidxs = floor((load_param.load.recs(2)-load_param.load.recs(1)+1) / param.get_heights.decimate_factor);
      rline0 = 1 + start_buffer;
      for adc_idx = 1:size(g_data,3)
        g_data(:,:,adc_idx) = fir_dec(g_data(:,:,adc_idx), param.get_heights.B_filter, ...
          param.get_heights.decimate_factor, rline0, Nidxs);
      end
      
      out_records.gps_time = fir_dec(out_records.gps_time, param.get_heights.B_filter, ...
        param.get_heights.decimate_factor, rline0, Nidxs);
      out_records.lat = fir_dec(out_records.lat, param.get_heights.B_filter, ...
        param.get_heights.decimate_factor, rline0, Nidxs);
      out_records.lon = fir_dec(out_records.lon, param.get_heights.B_filter, ...
        param.get_heights.decimate_factor, rline0, Nidxs);
      out_records.elev = fir_dec(out_records.elev, param.get_heights.B_filter, ...
        param.get_heights.decimate_factor, rline0, Nidxs);
      out_records.roll = fir_dec(out_records.roll, param.get_heights.B_filter, ...
        param.get_heights.decimate_factor, rline0, Nidxs);
      out_records.pitch = fir_dec(out_records.pitch, param.get_heights.B_filter, ...
        param.get_heights.decimate_factor, rline0, Nidxs);
      out_records.heading = fir_dec(out_records.heading, param.get_heights.B_filter, ...
        param.get_heights.decimate_factor, rline0, Nidxs);
    end
    
    gps_time = out_records.gps_time;
    lat = out_records.lat;
    lon = out_records.lon;
    elev = out_records.elev;
    roll = out_records.roll;
    pitch = out_records.pitch;
    heading = out_records.heading;
    
    %% 2. Extract surface values according to bin_rng
    layers(1).twtt = interp1(layers(1).gps_time, layers(1).twtt, gps_time(1,:));
    zero_bin = round(interp1(wfs(wf).time, 1:length(wfs(wf).time), layers(1).twtt,'linear','extrap'));
    start_bin = param.analysis.surf.bin_rng(1) + zero_bin;
    stop_bin = param.analysis.surf.bin_rng(end) + zero_bin;
    surf_vals = zeros(1+diff(param.analysis.surf.bin_rng([1 end])), size(g_data,2), size(g_data,3));
    surf_bins = start_bin;
    for rline = 1:size(g_data,2)
      start_bin0 = max(1,start_bin(rline));
      stop_bin0 = min(size(g_data,1),stop_bin(rline));
      out_bin0 = 1 + start_bin0-start_bin(rline);
      out_bin1 = size(surf_vals,1) - (stop_bin(rline)-stop_bin0);
      surf_vals(out_bin0:out_bin1,rline,:) = g_data(start_bin0:stop_bin0,rline,:);
    end

    %% 3. Save
    out_fn = fullfile(ct_filename_out(param, ...
      param.analysis.out_path, 'CSARP_noise'), ...
      sprintf('surf_img_%02d_%d_%d.mat',img_idx,param.load.recs(1),param.load.recs(end)));
    [out_fn_dir] = fileparts(out_fn);
    if ~exist(out_fn_dir,'dir')
      mkdir(out_fn_dir);
    end
    param_analysis = param;
    param_analysis.gps_source = records.gps_source;
    fprintf('  Saving outputs %s\n', out_fn);
    save(out_fn, 'surf_vals','surf_bins', 'wfs', 'gps_time', 'lat', ...
      'lon', 'elev', 'roll', 'pitch', 'heading', 'param_analysis', 'param_records');
  end
  
  if param.analysis.specular.en
    %% Compensate for elevation changes
    
    % Smooth elevation data since there seems to be small errors in it
    out_records.elev = sgolayfilt(out_records.elev, 3, 201, hanning(201));
    
    wf = abs(param.load.imgs{1}(1,1));
    adc = abs(param.load.imgs{1}(1,2));
    
    % Create frequency axis
    dt = wfs(wf).time{adc}(2) - wfs(wf).time{adc}(1);
    Nt = length(wfs(wf).time{adc});
    T = Nt*dt;
    df = 1/T;
    freq = wfs(wf).fc + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
    if any(strcmpi(param.radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3'}))
      % HACK: Sign error in pulse compression causes fc to be negative
      if isfield(param.radar.wfs,'fc_sign') && param.radar.wfs.fc_sign < 0
        freq_hack = -wfs(wf).fc + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
      else
        freq_hack = wfs(wf).fc + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
      end
    end

    % Correct all the data to a constant elevation (no zero padding is
    % applied so wrap around could be an issue for DDC data)
    g_data = cell2mat(g_data);
    for rline = 1:size(g_data,2)
      elev_dt = (out_records.elev(rline) - out_records.elev(1)) / (c/2);
      if any(strcmpi(param.radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3'}))
        g_data(:,rline) = ifft(ifftshift(fft(g_data(:,rline)),1) .* exp(j*2*pi*freq_hack*elev_dt));
      else
        g_data(:,rline) = ifft(fft(g_data(:,rline)) .* exp(j*2*pi*freq*elev_dt));
      end
    end
    
    %% Estimate slow-time coherence for every N range lines
    
    % Grab the peak values
    if ~isfield(param.analysis.specular,'min_bin') || isempty(param.analysis.specular.min_bin)
      param.analysis.specular.min_bin = 0e-6;
    end
    min_bin_idxs = find(wfs(wf).time{1} >= param.analysis.specular.min_bin,1);
    [max_value,max_idx_unfilt] = max(g_data(min_bin_idxs:end,:));
    max_idx_unfilt = max_idx_unfilt + min_bin_idxs(1) - 1;

    % Perform STFT (short time Fourier transform) (i.e. overlapping short FFTs in slow-time)
    H = spectrogram(double(max_value),hanning(param.analysis.specular.ave),param.analysis.specular.ave/2,param.analysis.specular.ave);

    % Since there may be a little slope in the ice, we sum the powers from
    % the lower frequency doppler bins rather than just taking DC. It seems to help
    % a lot to normalize by the sum of the middle/high-frequency Doppler bins.   A coherent/specular
    % surface will have high power in the low bins and low power in the high bins
    % so this ratio makes sense.
    peakiness = lp(max(abs(H(param.analysis.specular.signal_doppler_bins,:)).^2) ./ mean(abs(H(param.analysis.specular.noise_doppler_bins,:)).^2));
    
    if 0
      figure(1); clf;
      imagesc(lp(g_data))
      figure(2); clf;
      plot(peakiness)
      keyboard
    end
    
    % Threshold to find high peakiness range lines. (Note these are not
    % actual range line numbers, but rather indices into the STFT groups
    % of range lines.)
    good_rlines_mask = peakiness > param.analysis.specular.threshold;
    good_rlines = find(peakiness > param.analysis.specular.threshold);
    
    % Force there to be two good STFT groups in a row before storing
    % it to the specular file for deconvolution.
    good_rlines_idxs = diff(good_rlines) == 1;
    final_good_rlines = good_rlines(good_rlines_idxs);

    % Prepare outputs for file
    peakiness_rlines = round((1:length(peakiness)+0.5)*param.analysis.specular.ave/2);
    gps_time = records.gps_time(peakiness_rlines);
    lat = records.lat(peakiness_rlines);
    lon = records.lon(peakiness_rlines);
    elev = records.elev(peakiness_rlines);
    roll = records.roll(peakiness_rlines);
    pitch = records.pitch(peakiness_rlines);
    heading = records.heading(peakiness_rlines);

    deconv_forced = zeros(size(final_good_rlines));
    if isfield(param.analysis.specular,'gps_times') && ~isempty(param.analysis.specular.gps_times)
      %% Check to see if any of the forced GPS times are in this block
      for idx = 1:length(param.analysis.specular.gps_times)
        force_gps_time = param.analysis.specular.gps_times(idx);
        if records.gps_time(1) <= force_gps_time && records.gps_time(end) >= force_gps_time
          % This forced GPS time is in the block, find the peakiness block
          % closest to this time and force it to be included in final_good_rlines
          % if it is not already.
          [~,force_final_good_rline] = min(abs(gps_time - force_gps_time));
          match_idx = find(final_good_rlines == force_final_good_rline);
          if isempty(match_idx)
            final_good_rlines = [final_good_rlines force_final_good_rline];
            [final_good_rlines,new_idxs] = sort(final_good_rlines);
            deconv_forced(new_idxs(end)) = 1;
          else
            deconv_forced(match_idx) = 1;
          end
        end
      end
    end

    % Look through each good STFT group
    deconv_gps_time = [];
    deconv_mean = {};
    deconv_std = {};
    deconv_sample = {};
    deconv_twtt = [];
    deconv_DDC_Mt = [];
    for good_rline_idx = 1:length(final_good_rlines)
      % Get the specific STFT group we will be extracting an answer from
      final_good_rline = final_good_rlines(good_rline_idx);
      
      % Determine the center range line that this STFT group corresponds to
      center_rline = (final_good_rline+0.5)*param.analysis.specular.ave/2;
      
      fprintf('    SPECULAR %d %s (%s)!\n', center_rline, ...
        datestr(epoch_to_datenum(records.gps_time(center_rline)),'YYYYmmDD HH:MM:SS.FFF'), ...
        datestr(now));
      
      % Find the max values and correponding indices for all the range lines
      % in this group. Since we over-interpolate by Mt and the memory
      % requirements may be prohibitive, we do this in a loop
      % Enforce the same DDC filter in this group. Skip groups that have DDC filter swiches. 
      STFT_rlines = -param.analysis.specular.ave/4 : param.analysis.specular.ave/4-1;         
      if any(diff(img_Mt{1}(center_rline + STFT_rlines)))
        fprintf('    Including different DDC filters, skipped.\n');
        continue
      end
      Mt = 100;
      max_value = zeros(size(STFT_rlines));
      max_idx_unfilt = zeros(size(STFT_rlines));
      for offset_idx = 1:length(STFT_rlines)
        offset = STFT_rlines(offset_idx);
        oversampled_rline = interpft(g_data(:,center_rline+offset),size(g_data,1)*Mt);
        [max_value(offset_idx),max_idx_unfilt(offset_idx)] ...
          = max(oversampled_rline(min_bin_idxs(1)*Mt:end));
        max_idx_unfilt(offset_idx) = max_idx_unfilt(offset_idx) + min_bin_idxs(1)*Mt - 1;
      end
      
      % Filter the max and phase vectors
      max_idx = sgolayfilt(max_idx_unfilt/100,3,51);
      phase_corr = sgolayfilt(double(unwrap(angle(max_value))),3,51);

      %% Compensate range lines for amplitude, phase, and delay variance
      % in the peak value
      
      % Apply true time delay shift to flatten surface
      comp_data = ifft(fft(g_data(:,center_rline+STFT_rlines)) .* exp(j*2*pi*freq*max_idx*dt) );
      % Apply amplitude correction
      comp_data = comp_data .* repmat(1./abs(max_value), [Nt 1]);
      % Apply phase correction (compensating for phase from time delay shift)
      comp_data = comp_data .* repmat(exp(-j*(phase_corr + 2*pi*wfs(1).fc*max_idx*dt)), [Nt 1]);

      deconv_gps_time(end+1) = records.gps_time(center_rline);
      deconv_mean{end+1} = mean(comp_data,2);
      deconv_std{end+1} = std(comp_data,[],2);
      deconv_sample{end+1} = g_data(:,center_rline+1+param.analysis.specular.ave/4);
      deconv_twtt(:,end+1) = wfs(1).time{adc}(round(mean(max_idx)));
      deconv_DDC_Mt(end+1) = img_Mt{1}(center_rline);
    end
        
    wfs(wf).freq = freq;
    
    out_fn = fullfile(ct_filename_out(param, ...
      param.analysis.out_path, 'CSARP_noise'), ...
      sprintf('specular_%d_%d.mat',param.load.recs(1),param.load.recs(end)));
    [out_fn_dir] = fileparts(out_fn);
    if ~exist(out_fn_dir,'dir')
      mkdir(out_fn_dir);
    end
    param_analysis = param;
    fprintf('  Saving outputs %s\n', out_fn);
    save(out_fn, 'deconv_gps_time', 'deconv_mean', 'deconv_std','deconv_sample','deconv_twtt',...
      'deconv_DDC_Mt','deconv_forced','peakiness', 'wfs', 'gps_time', 'lat', ...
      'lon', 'elev', 'roll', 'pitch', 'heading', 'param_analysis', 'param_records');
  end
  
  if param.analysis.coh_ave.en
    g_data = cell2mat(g_data);
    % Implement memory efficient fft operations
    doppler = zeros(1,size(g_data,2));
    for rbin=1:size(g_data,1)
      doppler = doppler + abs(fft(conj(g_data(mod(rbin-1,size(g_data,1))+1,:)))).^2;
    end
    doppler = doppler*size(g_data,1);
    for rline=1:size(g_data,2)
      g_data(:,rline) = fft(g_data(:,rline));
    end
    
    % Do averaging
    rline0_list = 1:param.analysis.coh_ave.block_ave:size(g_data,2);
    for rline0_idx = 1:length(rline0_list)
      rline0 = rline0_list(rline0_idx);
      rlines = rline0 + (0:min(param.analysis.coh_ave.block_ave-1,size(g_data,2)-rline0));
      
      if strcmp(param.radar_name,'kuband') ...
          && (strcmp(param.season_name,'2009_Antarctica_DC8') ...
          || strcmp(param.season_name,'2011_Greenland_P3'))
        %% HACK for time variant noise floor
        % Create frequency axis
        dt = wfs.time(2) - wfs.time(1);
        Nt = length(wfs.time);
        T = Nt*dt;
        df = 1/T;
        freq = fftshift(wfs.fc + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).');
        if 0
          % 2009_Antarctica_DC8
          data_mean_removed = g_data(:,rlines) - repmat(mean(g_data,2),[1 numel(rlines)]);
          noise_power = sgolayfilt(double(lp(mean(abs(data_mean_removed).^2,2))), 2, 1001);
          plot(noise_power);
          save(sprintf('~/%s_kuband_noise_power.mat', param.season_name), 'noise_power','freq');
        elseif 0
          % 2011_Greenland_P3, 20110316_01
          data_mean_removed = g_data(:,rlines) - repmat(mean(g_data,2),[1 numel(rlines)]);
          noise_power = sgolayfilt(double(lp(mean(abs(data_mean_removed).^2,2))), 2, 501);
          noise_power(3200:4400) = NaN;
          plot(noise_power);
          save(sprintf('~/%s_kuband_noise_power.mat', param.season_name), 'noise_power','freq');
        end
        coh = load(sprintf('~/%s_kuband_noise_power.mat', param.season_name), 'noise_power','freq');
        coh.noise_power = interp1(coh.freq(~isnan(coh.noise_power)), ...
          coh.noise_power(~isnan(coh.noise_power)),freq,'nearest','extrap');
        good_samples = lp(g_data(:,rlines) - repmat(mean(g_data,2),[1 numel(rlines)])) ...
          < repmat(coh.noise_power+param.analysis.coh_ave.power_threshold,[1 numel(rlines)]);
      else
        good_samples = lp(g_data(:,rlines) - repmat(mean(g_data,2),[1 numel(rlines)])) < param.analysis.coh_ave.power_threshold;
      end
      
      if 0
        % Debug Plots for determining coh_ave.power_threshold
        figure(1); clf;
        imagesc(lp(g_data(:,rlines)));
        a1 = gca;
        figure(2); clf;
        imagesc(good_samples);
        colormap(gray);
        title('Good sample mask (black is thresholded)');
        a2 = gca;
        figure(3); clf;
        imagesc(lp(g_data(:,rlines) - repmat(mean(g_data,2),[1 numel(rlines)])) );
        a3 = gca;
        linkaxes([a1 a2 a3], 'xy');
      end
      
      coh_ave_samples(:,rline0_idx) = sum(good_samples,2);
      coh_ave(:,rline0_idx) = sum(g_data(:,rlines) .* good_samples,2) ./ coh_ave_samples(:,rline0_idx);
      gps_time(rline0_idx) = mean(records.gps_time(rlines));
      lat(rline0_idx) = mean(records.lat(rlines));
      lon(rline0_idx) = mean(records.lon(rlines));
      elev(rline0_idx) = mean(records.elev(rlines));
      roll(rline0_idx) = mean(records.roll(rlines));
      pitch(rline0_idx) = mean(records.pitch(rlines));
      heading(rline0_idx) = mean(records.heading(rlines));
    end
    
    time = wfs(wf).time;
    
    out_fn = fullfile(ct_filename_out(param, ...
      param.analysis.out_path, 'CSARP_noise'), ...
      sprintf('coh_noise_%d_%d.mat',param.load.recs(1),param.load.recs(end)));
    [out_fn_dir] = fileparts(out_fn);
    if ~exist(out_fn_dir,'dir')
      mkdir(out_fn_dir);
    end
    param_analysis = param;
    fprintf('  Saving outputs %s\n', out_fn);
    save(out_fn, 'coh_ave', 'coh_ave_samples', 'doppler', 'time', 'gps_time', 'lat', ...
      'lon', 'elev', 'roll', 'pitch', 'heading', 'param_analysis', 'param_records');
  end
  
end

success = true;

return;
