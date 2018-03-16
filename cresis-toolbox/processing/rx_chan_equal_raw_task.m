function [success] = rx_chan_equal_raw_task(param)
% [success] = rx_chan_equal_raw_task(param)
%
% Cluster task for rx_chan_equal_raw. Does the actual data loading
% and layer tracking.
%
% param = struct controlling the loading, processing, layer tracking,
%   and receiver equalization
%
% success = boolean which is true when the function executes properly
%   if a task fails before it can return success, then success will be
%   empty
%
% Author: John Paden
%
% See also rx_chan_equal_raw.m

global g_data;

physical_constants;
surfTimes = [];

records_fn = ct_filename_support(param,'','records');

if ~isfield(param.equal,'elev_correction') || isempty(param.equal.elev_correction)
  param.equal.elev_correction = false;
end

if ~param.equal.elev_correction ...
    && length(param.equal.B_filter) == param.equal.decimate_factor ...
    && all(param.equal.B_filter == 1)
  simple_firdec = true;
else
  simple_firdec = false;
end

if ~isfield(param.equal,'trim_vals') || isempty(param.equal.trim_vals)
  param.equal.trim_vals = [0 0];
end

if ~isfield(param.equal,'coh_noise_method') || isempty(param.equal.coh_noise_method)
  param.equal.coh_noise_method = 0;
end

if ~isfield(param.equal,'pulse_rfi')
  param.equal.pulse_rfi.en = 0;
end

if ~isfield(param.equal,'ft_dec')
  param.equal.ft_dec = 1;
end

if ~isfield(param.records,'file_version')
  param.records.file_version = [];
end

if ~isfield(param.equal,'layer_src') || isempty(param.equal.layer_src)
  % {'layerData',PATH,1}, {'records'}, {'ops',layer_name_string}, or numeric twtt
  param.equal.layer_src = {'layerData','layerData',1};
end

%% =====================================================================
% Determine which records to load with load_mcords_data
%
% Load records on either side of the current block, note if at the
% beginning or end of the segment.  Load with minimal presumming.

if simple_firdec
  load_param.load.recs(1) = param.load.recs(1);
  load_param.load.recs(2) = param.load.recs(2);
  records = read_records_aux_files(records_fn,load_param.load.recs);
  old_param_records = records.param_records;
else
  if mod(length(param.equal.B_filter)-1,2)
    error('Filter order must be even (e.g. fir1(EVEN_NUMBER,cutoff))');
  end
  filter_order = length(param.equal.B_filter) - 1;
  start_buffer = min(filter_order/2,param.load.recs(1)-1);
  load_param.load.recs(1) = param.load.recs(1)-start_buffer;
  load_param.load.recs(2) = param.load.recs(2)+filter_order/2;
  records = read_records_aux_files(records_fn,load_param.load.recs);
  old_param_records = records.param_records;
  stop_buffer = filter_order/2 - ((load_param.load.recs(2)-load_param.load.recs(1)+1) ...
    - length(records.lat));
end
old_param_records.gps_source = records.gps_source;

%% =====================================================================
% Collect waveform information into one structure
%  (used by load_RADARNAME_data)

if simple_firdec
  param.equal.presums = param.equal.decimate_factor;
else
  param.equal.presums = 1;
end
param.equal.pulse_comp = 1;
param.equal.combine_rx = 0;
if strcmpi(param.radar_name,'mcrds')
  [wfs,rec_data_size] = load_mcrds_wfs(records.settings, param, ...
    1:max(old_param_records.records.file.adcs), param.equal);
  load_param.load.rec_data_size = rec_data_size;
elseif any(strcmpi(param.radar_name,{'mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
  [wfs,rec_data_size] = load_mcords_wfs(records.settings, param, ...
    1:max(old_param_records.records.file.adcs), param.equal);
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
    1:max(old_param_records.records.file.adcs), param.equal);
end
% load_wf_cmd=['load_',param.radar_name,'_wfs(records.wfs,param,1:max(old_param_records.file.adcs),param.equal)'];
% [wfs,rec_data_size] = eval(load_wf_cmd);
load_param.wfs = wfs;

%% =====================================================================
% Load layer information

if iscell(param.equal.layer_src) && strcmpi(param.equal.layer_src{1},'records')
  % Don't need to do anything
  records.layer = records.surface;

elseif iscell(param.equal.layer_src) && strcmpi(param.equal.layer_src{1},'ops')
  % Load layer data from OPS
  
elseif iscell(param.equal.layer_src) && strcmpi(param.equal.layer_src{1},'layerData')
  layer_fn = param.equal.layer_src{2};
  layer_idx = param.equal.layer_src{3};
  
  %% Get the generic layer data path
  layer_path = fullfile(ct_filename_out(param,layer_fn,'',0));
  
  %% Load the current frame
  layer_fn = fullfile(layer_path,sprintf('Data_%s_%03d.mat',param.day_seg,param.proc.frm));
  layer = load(layer_fn);
  new_layer_gps_time = layer.GPS_time;
  layer = layer.layerData{layer_idx}.value{2}.data;
  
  %% Get the previous frame if necessary
  if records.gps_time(1) < new_layer_gps_time(1)-1
    layer_fn = fullfile(layer_path,sprintf('Data_%s_%03d.mat',param.day_seg,param.proc.frm-1));
    if exist(layer_fn,'file')
      tmp = load(layer_fn);
      new_layer_gps_time = [tmp.GPS_time new_layer_gps_time];
      layer = [tmp.layerData{layer_idx}.value{2}.data layer];
    end
  end
  %% Get the next frame if necessary
  if records.gps_time(end) > new_layer_gps_time(end)+1
    layer_fn = fullfile(layer_path,sprintf('Data_%s_%03d.mat',param.day_seg,param.proc.frm+1));
    if exist(layer_fn,'file')
      tmp = load(layer_fn);
      new_layer_gps_time = [new_layer_gps_time tmp.GPS_time];
      layer = [layer tmp.layerData{layer_idx}.value{2}.data];
    end
  end
  
  %% Since layer files may have overlapping data, sort it
  [new_layer_gps_time new_layer_idxs] = sort(new_layer_gps_time);
  layer = layer(new_layer_idxs);
  
  % Remove NaN
  good_mask = ~isnan(new_layer_gps_time);
  new_layer_gps_time = new_layer_gps_time(good_mask);
  layer = layer(good_mask);
  
  %% Do the interpolation from layer data file to records gps time
  layer = interp1(new_layer_gps_time,layer,records.gps_time,'linear','extrap');
  
  %% Fill in NaN gaps
  records.layer_mask = isfinite(layer);
  records.layer = interp_finite(layer,1);
end

%% =====================================================================
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
if any(strcmpi(param.radar_name,{'mcords','mcords2','mcords3','mcords4','mcords5','accum2'}))
  boards = adc_to_board(param.radar_name,records.param_records.records.file.adcs);
  for idx = 1:length(param.load.adcs)
    adc = param.load.adcs(idx);
    adc_idx = find(old_param_records.records.file.adcs == param.load.adcs(idx));
    if isempty(adc_idx)
      error('ADC %d not present in records file\n', param.load.adcs(idx));
    end
    
    % Just get the file-information for the records we need
    board = adc_to_board(param.radar_name,adc);
    board_idx = find(board == boards);
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
  load_param.load.offset = records.offset;
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

load_param.proc.trim_vals           = param.equal.trim_vals;
load_param.proc.pulse_comp          = param.equal.pulse_comp;
load_param.proc.ft_dec              = param.equal.ft_dec;
load_param.proc.ft_wind             = param.equal.ft_wind;
load_param.proc.ft_wind_time        = param.equal.ft_wind_time;
load_param.proc.presums             = param.equal.presums;
load_param.proc.combine_rx          = param.equal.combine_rx;
load_param.proc.pulse_rfi           = param.equal.pulse_rfi;
load_param.proc.coh_noise_method    = param.equal.coh_noise_method;
load_param.proc.coh_noise_arg       = param.equal.coh_noise_arg;

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

for img = 1:length(param.load.imgs)
  
  if isnumeric(param.equal.layer_src)
    wf = param.equal.imgs{1}(1);
    records.layer = param.equal.layer_src * ones(size(records.gps_time));
    records.layer_mask = ones(size(records.gps_time));
  end

  out_records = records;

  %% Load data into g_data using load_mcords_data
  load_param.load.imgs = param.load.imgs(img);
  % Determine combination times when multiple wf-adc pairs are being loaded
  % to form a single range line
  if size(load_param.load.imgs{1},2) == 2
    load_param.load.wf_adc_comb.en = 0;
  else
    load_param.load.wf_adc_comb.en = 1;
    wf1 = load_param.load.imgs{1}(1,1);
    wf2 = load_param.load.imgs{1}(1,3);
    % t1 = time to switch from wf1 to wf2
    %wf_adc_layer = fir_dec(records.layer, param.equal.decimate_factor);
    wf_adc_layer = records.layer;
    t1 = eval(param.equal.qlook.wf_adc_comb{1});
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
    wfs.time = load_fmcw_data(load_param,out_records.gps_time);
  end
  
  %% Remove overlap data
  recs_keep = 1+param.load.recs_keep(1)-load_param.load.recs(1) ...
    : length(out_records.lat)+param.load.recs_keep(end)-load_param.load.recs(end);
  out_records.lat = out_records.lat(recs_keep);
  out_records.lon = out_records.lon(recs_keep);
  out_records.elev = out_records.elev(recs_keep);
  out_records.roll = out_records.roll(recs_keep);
  out_records.pitch = out_records.pitch(recs_keep);
  out_records.heading = out_records.heading(recs_keep);
  out_records.gps_time = out_records.gps_time(recs_keep);
  out_records.layer = out_records.layer(recs_keep);
  out_records.layer_mask = out_records.layer_mask(recs_keep);
  if simple_firdec
    recs_keep = 1+(param.load.recs_keep(1)-load_param.load.recs(1))/param.equal.decimate_factor ...
      : size(g_data,2)+(param.load.recs_keep(end)-load_param.load.recs(end))/param.equal.decimate_factor;
  end
  g_data = g_data(:,recs_keep,:);
  
  %% Remove coherent noise
  if param.equal.coh_noise_method == 1 && ~any(strcmpi(param.radar_name,{'kuband','snow','kuband2','snow2','kuband3','snow3'}))
    for wf_adc_idx = 1:size(g_data,3)
      g_data(:,:,wf_adc_idx) = g_data(:,:,wf_adc_idx) - repmat( mean(g_data(:,:,wf_adc_idx),2) , ...
        [1 size(g_data,2)]);
    end
  end
  
  %% FIR Decimate
  if simple_firdec
%     if img == 1
      out_records.gps_time = fir_dec(out_records.gps_time, param.equal.decimate_factor);
      out_records.lat = fir_dec(out_records.lat, param.equal.decimate_factor);
      out_records.lon = fir_dec(out_records.lon, param.equal.decimate_factor);
      out_records.elev = fir_dec(out_records.elev, param.equal.decimate_factor);
      out_records.roll = fir_dec(out_records.roll, param.equal.decimate_factor);
      out_records.pitch = fir_dec(out_records.pitch, param.equal.decimate_factor);
      out_records.heading = fir_dec(out_records.heading, param.equal.decimate_factor);
      out_records.layer = fir_dec(out_records.layer, param.equal.decimate_factor);
      out_records.layer_mask = fir_dec(out_records.layer_mask, param.equal.decimate_factor);
%     end
    
  else
    % Check for edge conditions that caused not enough data to be loaded
    % in the case where a segment starts and stops.
    % If not enough data was loaded, modify the filter coefficient
    % normalization so that it handles this
    Nidxs = floor((load_param.load.recs(2)-load_param.load.recs(1)+1) / param.equal.decimate_factor);
    rline0 = 1 + start_buffer;
    g_data = fir_dec(g_data, param.equal.B_filter, ...
      param.equal.decimate_factor, rline0, Nidxs);

%     if img == 1
      out_records.gps_time = fir_dec(out_records.gps_time, param.equal.B_filter, ...
        param.equal.decimate_factor, rline0, Nidxs);
      out_records.lat = fir_dec(out_records.lat, param.equal.B_filter, ...
        param.equal.decimate_factor, rline0, Nidxs);
      out_records.lon = fir_dec(out_records.lon, param.equal.B_filter, ...
        param.equal.decimate_factor, rline0, Nidxs);
      out_records.elev = fir_dec(out_records.elev, param.equal.B_filter, ...
        param.equal.decimate_factor, rline0, Nidxs);
      out_records.roll = fir_dec(out_records.roll, param.equal.B_filter, ...
        param.equal.decimate_factor, rline0, Nidxs);
      out_records.pitch = fir_dec(out_records.pitch, param.equal.B_filter, ...
        param.equal.decimate_factor, rline0, Nidxs);
      out_records.heading = fir_dec(out_records.heading, param.equal.B_filter, ...
        param.equal.decimate_factor, rline0, Nidxs);
      out_records.layer = fir_dec(out_records.layer, param.equal.B_filter, ...
        param.equal.decimate_factor, rline0, Nidxs);
      out_records.layer_mask = fir_dec(out_records.layer_mask, param.equal.B_filter, ...
        param.equal.decimate_factor, rline0, Nidxs);
%     end
  end
  
  wf = abs(param.load.imgs{img}(1,1));
  Time = wfs(wf).time;
  
  %% Motion compensation
  % ======================================================================
  mocomp_param.type = param.equal.mocomp_type;
  mocomp_param.tx_weights = param.radar.wfs(wf).tx_weights;
  mocomp_param.season_name = param.season_name;
  mocomp_param.radar_name = param.radar_name;
  mocomp_param.gps_source = out_records.gps_source;
  for wf_adc_idx = 1:size(g_data,3)
    wf = abs(param.load.imgs{img}(wf_adc_idx,1));
    adc = abs(param.load.imgs{img}(wf_adc_idx,2));
    mocomp_param.rx = param.radar.wfs(wf).rx_paths(adc);
    
    % drange = change in range (positive is longer range)
    drange = basic_motion_comp(mocomp_param,param.equal.lever_arm_fh,out_records.roll, ...
      out_records.pitch,out_records.heading,out_records.lat,out_records.lon,out_records.elev);
    % dtime = Convert to time (in air), positive is longer range/time-delay
    dtime = 2*drange/3e8;
    Nt = size(g_data,1);
    Nx = size(g_data,2);
    if 0
      figure(1); clf;
      plot(dtime*1e9);
      title('adc %d');
      pause;
    end
    % Time shift data in the frequency domain
    %   A positive dtime implies a larger negative phase delay (since
    %   delay is negative/lagging phase)
    g_data(:,:,wf_adc_idx) = ifft(fft(g_data(:,:,wf_adc_idx)) ...
      .*exp(-1i*2*pi*repmat(wfs(wf).freq,1,Nx).*repmat(dtime,Nt,1)));
  end
  
  %% Convert layer twtt to range bins
  % =======================================================================
  % retrack_bins is usually -N:N vector (0 causes no retracking to be done)
  %   It specifies the range of bins relative to the current bin that will
  %   be examined during retracking.
  retrack_bins = param.equal.layer_retrack;
  
  layer_bins = interp1(Time,1:length(Time),out_records.layer);
  % Insure that retracker won't have indexing problems
  layer_bins(out_records.layer < Time( 1-retrack_bins(1) )) = 1-retrack_bins(1);
  out_records.layer_mask(out_records.layer < Time( 1-retrack_bins(1) )) = 0;
  layer_bins(out_records.layer > Time(end-retrack_bins(end))) = length(Time)-retrack_bins(end);
  out_records.layer_mask(out_records.layer > Time(end-retrack_bins(end))) = 0;
  % Round since we use this as the index
  layer_bins = round(layer_bins);

  %% Retrack layer
  % =======================================================================
  data_for_tracking = fir_dec(abs(mean(g_data,3)).^2,ones(1,9),1);
  for rline=1:size(data_for_tracking,2)
    [layer_bins_vals(rline),offset] = max(data_for_tracking(layer_bins(rline)+retrack_bins,rline));
    layer_bins(rline) = layer_bins(rline) + retrack_bins(1) - 1 + offset;
  end
  
  if 0
    for chan = 1:size(g_data,3)
      figure(1); clf;
      imagesc(angle(g_data(:,:,chan) .* conj(g_data(:,:,12))))
      hold on
      plot(layer_bins)
      hold off
      ylim([150 200]);
      
      figure(2); clf;
      imagesc(lp(g_data(:,:,chan)))
      hold on
      plot(layer_bins)
      hold off
      ylim([150 200]);
      
      pause;
    end
  end
  
  %% Cross correlation to determine recommended time, phase, and amplitude offsets
  % =======================================================================
  ref_idx = param.equal.ref_wf_adc_idx;
  
  ref_bins = param.equal.ref_bins(1):param.equal.ref_bins(2);
  search_bins = param.equal.search_bins(1)+param.equal.ref_bins(1) : param.equal.search_bins(2)+param.equal.ref_bins(2);
  zero_padding_offset = length(search_bins) - length(ref_bins);
  Hcorr_wind = hanning(length(ref_bins));
  clear peak_val peak_offset;
  
  rlines = 1:size(g_data,2); % We always use all the range lines... if not, this would be the variable to set
  peak_val = zeros(size(g_data,3),length(rlines));
  peak_val_corr = zeros(size(g_data,3),length(rlines));
  peak_offset = zeros(size(g_data,3),length(rlines));
  
  for adc_idx = 1:size(g_data,3)
    for rline_idx = 1:size(g_data,2)
      rline = rlines(rline_idx);
      
      bin = layer_bins(rline_idx);
      
      if bin+search_bins(1) < 1
        search_bins_trunc = 1-bin : search_bins(end);
      elseif  bin+search_bins(end) > Nt
        search_bins_trunc = search_bins(1) : Nt-bin;
      else
        search_bins_trunc = search_bins;
      end
       
      if bin+ref_bins(1) < 1
        ref_bins_trunc = 1-bin : ref_bins(end);
        Hcorr_wind_trunc = Hcorr_wind(2-ref_bins(1)-bin:end);
      elseif  bin+ref_bins(end) > Nt
        ref_bins_trunc = ref_bins(1) : Nt-bin;
        Hcorr_wind_trunc = Hcorr_wind(1:end-ref_bins(end) + (Nt-bin));
      else
        ref_bins_trunc = ref_bins;
        Hcorr_wind_trunc = Hcorr_wind;
      end
      
      [corr_out,lags] = xcorr(g_data(layer_bins(rline_idx)+search_bins_trunc,rline,adc_idx), ...
        g_data(layer_bins(rline_idx)+ref_bins_trunc,rline,ref_idx) .* Hcorr_wind_trunc);
      corr_int = interpft(corr_out,param.equal.Mt*length(corr_out));
      [peak_val_corr(adc_idx,rline_idx) peak_offset(adc_idx,rline_idx)] = max(corr_int);
      peak_val(adc_idx,rline_idx) = g_data(layer_bins(rline_idx),rline,adc_idx) ...
        .* conj(g_data(layer_bins(rline_idx),rline,ref_idx));
      peak_offset(adc_idx,rline_idx) = (peak_offset(adc_idx,rline_idx)-1)/param.equal.Mt+1 ...
        + ref_bins_trunc(1) + search_bins_trunc(1) - 1 - zero_padding_offset;
    end
  end
  
  if 0
    figure(1); clf;
    h = plot(angle(peak_val ./ repmat(peak_val(ref_idx,:),[size(peak_val,1) 1])).'*180/pi,'.');
    legend(h,{'1','2','3','4','5','6'});
    grid on;
    
    hold on;
    h = plot(angle(peak_val_corr ./ repmat(peak_val_corr(ref_idx,:),[size(peak_val,1) 1])).'*180/pi);
    hold off;
    set(h,'LineStyle','--');
    legend(h,{'1','2','3','4','5','6'});
    grid on;
  end
  
  peak_offset = peak_offset - repmat(peak_offset(ref_idx,:),[size(peak_offset,1),1]);
  dt = (wfs(wf).time(2)-wfs(wf).time(1));
  
  %% Roll Estimation
  % =======================================================================
  if isempty(param.equal.lever_arm_fh)
    roll_est_theta = [];
    roll_est_val = [];
  else
    param.roll_est.bin_rng = 0;
    param.roll_est.rline_rng = -5:5;
    param.roll_est.Nsig = 1;
    y_offset = zeros(size(g_data,3),1);
    z_offset = zeros(size(g_data,3),1);
    lever_arm_param.season_name = param.season_name;
    lever_arm_param.radar_name = ct_output_dir(param.radar_name);
    lever_arm_param.gps_source = out_records.gps_source;
    for wf_adc_idx = 1:size(g_data,3)
      wf = abs(param.equal.imgs{img}(wf_adc_idx,1));
      adc = param.equal.imgs{img}(wf_adc_idx,2);
      mocomp_param.rx = param.radar.wfs(wf).rx_paths(adc);
      phase_center = param.equal.lever_arm_fh(lever_arm_param, param.radar.wfs(wf).tx_weights, mocomp_param.rx);
      y_offset(wf_adc_idx) = -phase_center(2);
      z_offset(wf_adc_idx) = -phase_center(3);
    end
    
    [theta,sv] = array_proc_sv(256,wfs(wf).freq(1), y_offset, z_offset);
    theta = fftshift(theta);
    sv = fftshift(sv,2);
    
    Nt = size(g_data,1);
    Nx = size(g_data,2);
    Nc = size(g_data,3);
    for rline = 1:length(rlines)
      bin = layer_bins(rline);
      
      if rline+param.roll_est.rline_rng(1) < 1
        rline_rng = 1-rline : param.roll_est.rline_rng(end);
      elseif  rline+param.roll_est.rline_rng(end) > Nx
        rline_rng = param.roll_est.rline_rng(1) : Nx-rline;
      else
        rline_rng = param.roll_est.rline_rng;
      end
      
      if bin+param.roll_est.bin_rng(1) < 1
        bin_rng = 1-bin : param.roll_est.bin_rng(end);
      elseif  bin+param.roll_est.bin_rng(end) > Nt
        bin_rng = param.roll_est.bin_rng(1) : Nt-bin;
      else
        bin_rng = param.roll_est.bin_rng;
      end
      
      dataSample = g_data(bin+bin_rng,rline+rline_rng,:);
      dataSample = reshape(dataSample,[length(bin_rng)*length(rline_rng) Nc]).';
      
      Rxx = 1/size(dataSample,1) * (dataSample * dataSample');
      [V,D] = eig(Rxx);
      eigenVals = diag(D);
      [eigenVals noiseIdxs] = sort(eigenVals);
      noiseIdxs = noiseIdxs(1:end-param.roll_est.Nsig);
      music_pattern = mean(abs(sv'*V(:,noiseIdxs)).^2,2);
      [roll_est_val(rline),roll_est_idx] = min(music_pattern);
      roll_est_theta(rline) = theta(roll_est_idx);
    end
  end
  
  %% Save output
  % =======================================================================
  gps_time  =out_records.gps_time;
  gps_source =out_records.gps_source;
  roll = out_records.roll;
  layer_mask = out_records.layer_mask;
  
  fn = fullfile(ct_filename_out(param, ...
    param.equal.out_path, 'CSARP_equal'), ...
    sprintf('equal_%03d',param.proc.frm), ...
    sprintf('blk_%d_img_%02d.mat', param.proc.block, img));
  [fn_dir] = fileparts(fn);
  if ~exist(fn_dir,'dir')
    mkdir(fn_dir);
  end
  param_records = old_param_records;
  param_equal = param;
  fprintf('  Saving %s\n', fn);
  save(fn, 'layer_bins', 'layer_mask', 'layer_bins_vals','peak_val_corr', ...
    'peak_offset','peak_val','roll_est_theta','roll_est_val','gps_time', ...
    'gps_source','roll','param_equal', 'param_records','wfs');
  
end

success = true;

return;
