function [success] = qlook_task(param)
% [success] = qlook_task(param)
%
% Cluster task for quick look processing. Does the actual data loading
% and surface tracking.
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
% .qlook = structure controlling qlook processing
%  .radar_name = name of radar string
%  .season_name = name of mission string
%  .day_seg = day-segment string
%  
%  qlook fields used by load_mcords_data.m (see that function for details)
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
%  qlook fields for post processing
%  .roll_correction = boolean, whether or not to apply roll phase correction
%  .lever_arm_fh = string containing function name
%  .elev_correction = boolean, whether or not to apply elevation phase correction
%  .B_filter = double vector, FIR filter coefficients to apply before
%    decimating, this function loads data before and after this frame
%    (if available) to avoid transients at the beginning and end
%  .decimate_factor = positive integer, decimation rate
%  .inc_B_filter: double vector, FIR filter coefficients to apply before
%    incoherent average decimation. If not defined or empty, then
%    inc_B_filter is set to ones(1,inc_ave)/inc_ave.
%  .inc_ave = integer scalar, number of incoherent averages to apply
%    (also decimates by this number). If set to < 1, complex data are
%    returned.  Setting to 1 causes the data to be power detected (i.e.
%    become incoherent), but no averaging is done.
%
%  .surf = qlook structure controlling surface tracking
%   .en = boolean, whether or not to apply surface tracking
%   .wf_idx = positive integer, which waveform in the wfs list to apply surface tracking on
%   .min_bin = double scalar, the minimum range time that the surface can be tracked to.
%     This is used to keep the surface tracking routine from picking up the
%     feedthrough.  It requires a minimum elevation AGL.
%   .manual = boolean, whether or not to enable the manual tracking
%     interface.  Generally better to let the automated routine run, fix in
%     picker, and then update records (so surf.manual is mostly for debugging)
%
%  .qlook = qlook structure controlling quick look generation
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
%
% Author: John Paden
%
% See also qlook.m

global g_data;
g_data = [];

physical_constants;

records_fn = ct_filename_support(param,'','records');

[output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);

%% Load records
% =========================================================================

% 1. Determine which records to load, load records on either side of the
%    current block
buffer = floor((length(param.qlook.B_filter)-1)/2);
start_buffer = min(buffer,param.load.recs(1)-1); % Ensures that records before start are not requested
param.load.recs(1) = param.load.recs(1)-start_buffer;
param.load.recs(2) = param.load.recs(2)+buffer; % read_records will truncate if buffer extends beyond end
records = read_records_aux_files(records_fn,param.load.recs);
stop_buffer = buffer + length(records.gps_time) - (diff(param.load.recs)+1);

% 2. Store the parameters that were used to create the records file
old_param_records = records.param_records;
old_param_records.gps_source = records.gps_source;

%% Load layer data
% =========================================================================

if ~isempty(param.qlook.surf_layer)
  new_layer = opsLoadLayers(param,param.qlook.surf_layer);
end

if ~isempty(param.qlook.bottom_layer)
  new_layer = opsLoadLayers(param,param.qlook.bottom_layer);
end

save('/tmp/qlook_test.mat')
keyboard
load('/tmp/qlook_test.mat')

% =====================================================================
% Collect waveform information into one structure
[wfs,state] = data_load_wfs(param,records);

% =====================================================================
% Load data
[hdr,data] = data_load(param,records,wfs);

% =====================================================================
% Load and process each image separately
%
% For each image:
% 1. Load receiver data separately (minimal presumming)
% 2. Remove coherent noise (slow-time mean removal)
% 3. Apply roll correction
% 4. Combine receivers
% 5. Apply elevation compensation
% 6. FIR decimate the data
% =====================================================================

for img = 1:length(param.load.imgs)
  % Setup roll correction
  if param.qlook.roll_correction
    if isempty(param.qlook.lever_arm_fh)
      error('lever_arm_fh must be defined if roll_correction is enabled');
    end
    
    trajectory_param = struct('gps_source',records.gps_source, ...
      'season_name',param.season_name,'radar_name',param.radar_name,'rx_path', 0, ...
      'tx_weights', [], 'lever_arm_fh', param.qlook.lever_arm_fh);
    ref = trajectory_with_leverarm(records,trajectory_param);

    lever_arm_fh = param.qlook.lever_arm_fh;
    % Setup motion compensation (roll removal)
    radar_lever_arm = zeros(3,size(param.load.imgs{img},1));
    for wf_adc_idx = 1:size(param.load.imgs{img},1)
      wf = abs(param.load.imgs{img}(wf_adc_idx,1));
      adc = abs(param.load.imgs{img}(wf_adc_idx,2));
      radar_lever_arm(:,wf_adc_idx) = lever_arm_fh(trajectory_param,wfs(wf).tx_weights,wfs(wf).rx_paths(adc));
    end
  end
  
  % Default values to use
  wf = abs(param.load.imgs{img}(1,1));
  adc = abs(param.load.imgs{img}(1,2));
  lambda_fc = c/wfs(wf).fc;

  %% Compute trajectory using GPS/INS data and the lever arm
  % out_records = motion compensated data
  % records = original record information (not to be used below this section)
  if isempty(param.qlook.lever_arm_fh)
    out_records = records;
  else
    trajectory_param = struct('gps_source',old_param_records.gps_source, ...
          'season_name',param.season_name,'radar_name',param.radar_name, ...
          'rx_path', wfs(wf).rx_paths(adc), ...
      'tx_weights', wfs(wf).tx_weights, 'lever_arm_fh', param.qlook.lever_arm_fh);
    for tmp_wf_adc_idx = 2:size(param.load.imgs{1},1)
      tmp_wf = abs(param.load.imgs{1}(tmp_wf_adc_idx,1));
      tmp_adc = abs(param.load.imgs{1}(tmp_wf_adc_idx,2));
      trajectory_param.rx_path(tmp_wf_adc_idx) = wfs(tmp_wf).rx_paths(tmp_adc);
    end
    out_records = trajectory_with_leverarm(records,trajectory_param);
  end
  
  %% Load data into g_data using data_load.m
  
    [img_time,img_valid_rng,img_deconv_filter_idx,img_Mt] = load_fmcw_data(load_param,out_records);
    
  if ~exist('deconv_filter_idx','var')
    deconv_filter_idx = NaN(1,size(g_data,2));
  end
  
  if 0
    % Deconvolution Test Code
    keyboard
    dd = load('/mnt/products/ref_Tpd_3us_adc_1.mat');
    %                 ee = fft(dd.ref(1:20:end),length(g_data(:,1,1)));
    dd.ref = ifft(dd.fref);
    ee = interpft(fft(dd.ref(1:20:end)),length(g_data(:,80,1)));
    plot(lp(ee))
    ff = g_data(:,80,1);
    
    figure(1); clf;
    plot(lp(ff));
    hold on
    plot(lp(ifft(fft(ff) .* ee)),'r')
    hold off;
    
    ff = g_data(:,:,1);
    
    figure(1); clf;
    imagesc(lp(ff));
    figure(2); clf;
    imagesc(lp(ifft(fft(ff) .* repmat(ee,[1 size(ff,2)]))))
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
  out_records.surface = out_records.surface(recs_keep);
  if simple_firdec
    recs_keep = 1+(param.load.recs_keep(1)-load_param.load.recs(1))/param.qlook.decimate_factor ...
      : size(g_data,2)+(param.load.recs_keep(end)-load_param.load.recs(end))/param.qlook.decimate_factor;
  end
  g_data = g_data(:,recs_keep,:);
  deconv_filter_idx = deconv_filter_idx(:,recs_keep);
  
  %% Remove coherent noise
  if param.qlook.coh_noise_method && ~any(strcmpi(radar_name,{'kuband','kuband2','kuband3','kaband3','mcords5','snow','snow2','snow3','snow5','snow8'}))
    
    if param.qlook.coh_noise_method == 3 && isempty(param.qlook.coh_noise_arg)
      param.qlook.coh_noise_arg = 255;
    end
    
    % Remove only the DC Doppler component
    for wf_adc_idx = 1:size(g_data,3)
      if param.qlook.coh_noise_method == 1
        g_data(:,:,wf_adc_idx) = bsxfun(@minus, g_data(:,:,wf_adc_idx), ...
          mean(g_data(:,:,wf_adc_idx),2));
      elseif param.qlook.coh_noise_method == 3
        g_data(:,:,wf_adc_idx) = bsxfun(@minus, g_data(:,:,wf_adc_idx), ...
          fir_dec(g_data(:,:,wf_adc_idx),hanning(param.qlook.coh_noise_arg).'/(param.qlook.coh_noise_arg/2+0.5),1));
      end
    end
    
  end

  %% Roll compensation
  if param.qlook.roll_correction
    % Apply roll-only motion compensation
    for wf_adc_idx = 1:size(g_data,3)
      wf = abs(param.load.imgs{img}(wf_adc_idx,1));
      adc = abs(param.load.imgs{img}(wf_adc_idx,2));
      rx = wfs(wf).rx_paths(adc);
      for rline = 1:size(g_data,2)
        drange = radar_lever_arm(2,wf_adc_idx) * -tan(out_records.roll(rline));
        dphase = drange * 2 * 2*pi/lambda_fc;
        g_data(:,rline,wf_adc_idx) = g_data(:,rline,wf_adc_idx) * exp(1j*dphase);
      end
    end
    g_data = mean(g_data,3);
  end
  
  %% Elevation compensation
  if param.qlook.elev_correction && ~simple_firdec
    % Remove elevation variations (just a phase shift, not a time shift)
    %  - With simple_firdec there is no point in elevation compensation
    %    because the along-track averaging has already been done
    drange = out_records.elev-mean(out_records.elev);
    dphase = drange * 2 * 2*pi/lambda_fc;
    for rline = 1:size(g_data,2)
      g_data(:,rline) = g_data(:,rline) .* exp(1j*dphase(rline));
    end
  end

  %% FIR Decimate
  if simple_firdec
%     if img == 1
      out_records.gps_time = fir_dec(out_records.gps_time, param.qlook.decimate_factor);
      out_records.lat = fir_dec(out_records.lat, param.qlook.decimate_factor);
      out_records.lon = fir_dec(out_records.lon, param.qlook.decimate_factor);
      out_records.elev = fir_dec(out_records.elev, param.qlook.decimate_factor);
      out_records.roll = fir_dec(out_records.roll, param.qlook.decimate_factor);
      out_records.pitch = fir_dec(out_records.pitch, param.qlook.decimate_factor);
      out_records.heading = fir_dec(out_records.heading, param.qlook.decimate_factor);
%     end
    
  else
    % Check for edge conditions that caused not enough data to be loaded
    % in the case where a segment starts and stops.
    % If not enough data was loaded, modify the filter coefficient
    % normalization so that it handles this
    Nidxs = floor((param.load.recs(2)-param.load.recs(1)+1) / param.qlook.decimate_factor);
    rline0 = 1 + start_buffer;
    g_data = fir_dec(g_data, param.qlook.B_filter, ...
        param.qlook.decimate_factor, rline0, Nidxs);

%     if img == 1
      out_records.gps_time = fir_dec(out_records.gps_time, param.qlook.B_filter, ...
        param.qlook.decimate_factor, rline0, Nidxs);
      out_records.lat = fir_dec(out_records.lat, param.qlook.B_filter, ...
        param.qlook.decimate_factor, rline0, Nidxs);
      out_records.lon = fir_dec(out_records.lon, param.qlook.B_filter, ...
        param.qlook.decimate_factor, rline0, Nidxs);
      out_records.elev = fir_dec(out_records.elev, param.qlook.B_filter, ...
        param.qlook.decimate_factor, rline0, Nidxs);
      out_records.roll = fir_dec(out_records.roll, param.qlook.B_filter, ...
        param.qlook.decimate_factor, rline0, Nidxs);
      out_records.pitch = fir_dec(out_records.pitch, param.qlook.B_filter, ...
        param.qlook.decimate_factor, rline0, Nidxs);
      out_records.heading = fir_dec(out_records.heading, param.qlook.B_filter, ...
        param.qlook.decimate_factor, rline0, Nidxs);
%     end
  end
  
  if 0
    % Enable this if-statement only for debugging
    figure(img); clf;
    imagesc([1 size(g_data,2)],wfs(img).time, ...
      lp(g_data));
    keyboard
  end
  
  %% Apply incoherent averaging with decimation
  if param.qlook.inc_ave >= 1
    data_incoh = [];
    for adc_idx = 1:size(g_data,3)
      data_incoh(:,:,adc_idx) = fir_dec(fir_dec(abs(g_data(:,:,adc_idx)).^2,param.qlook.inc_B_filter,1), param.qlook.inc_ave);
    end
  end
  
  if 0
    % Undo elevation phase correction
    for rline = 1:size(g_data,2)
      drange = out_records.elev-mean(out_records.elev);
      dphase = drange * 2 * 2*pi/lambda_fc;
      g_data(:,rline) = g_data(:,rline) .* exp(-1j*dphase(rline));
    end
    
    % Leading edge detector
    data_lp = lp(data_incoh);
    [maxVal maxIdx] = max(data_lp);
    sizeBins = zeros(size(maxVal));
    for rline = 1:size(data_lp,2)
      surfBins(rline) = find(data_lp(:,rline) > maxVal(rline)-10,1);
      [tmp maxOffset] = max(data_lp(surfBins(rline)+(0:4),rline));
      surfBins2(rline) = surfBins(rline) + maxOffset-1;
    end
    surfBins3 = medfilt1(surfBins2,5);
    surfBins4 = reshape(repmat(surfBins3,[param.qlook.inc_ave 1]),[1 size(g_data,2)]);
    
    % Phase and elevation comparison (used for finding UTC time offset
    % errors)
    physical_constants;
    
    figure(1); clf;
    plot(out_records.elev);
    hold on;
    plot(-unwrap(angle(g_data(179,:))) * lambda_fc/(2*pi)/2 + 2463-0.1,'r');
    hold off;
    keyboard
  end
  
  % Save quick look results
  
  if param.qlook.inc_ave >= 1
    Data = data_incoh;
  else
    Data = g_data;
  end
  
  Time = wfs(wf).time;
  
  GPS_time = fir_dec(out_records.gps_time,param.qlook.inc_ave);
  Latitude = fir_dec(out_records.lat,param.qlook.inc_ave);
  Longitude = fir_dec(out_records.lon,param.qlook.inc_ave);
  Elevation = fir_dec(out_records.elev,param.qlook.inc_ave);
  Roll = fir_dec(out_records.roll,param.qlook.inc_ave);
  Pitch = fir_dec(out_records.pitch,param.qlook.inc_ave);
  Heading = fir_dec(out_records.heading,param.qlook.inc_ave);
  deconv_filter_idx = fir_dec(deconv_filter_idx,param.qlook.inc_ave);
  
  out_fn_name = sprintf('qlook_img_%02d_%d_%d.mat',img,param.load.recs_keep(1),param.load.recs_keep(end));
  out_fn_dir = fullfile(ct_filename_out(param, ...
    param.qlook.out_path, 'CSARP_qlook'), ...
    sprintf('ql_data_%03d_01_01',param.load.frm));
  out_fn = fullfile(out_fn_dir,out_fn_name);
  fprintf('  Save %s\n', out_fn);
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  param_records = old_param_records;
  param_qlook = param;
  custom = [];
  if any(~isnan(deconv_filter_idx))
    custom.deconv_filter_idx = deconv_filter_idx;
  end
  clear deconv_filter_idx;
  save(out_fn,'-v7.3', 'Data', 'Time', 'GPS_time', 'Latitude', ...
    'Longitude', 'Elevation', 'Roll', 'Pitch', 'Heading', 'param_qlook', 'param_records','custom');
  
end

fprintf('%s done %s\n', mfilename, datestr(now));

success = true;

return;
