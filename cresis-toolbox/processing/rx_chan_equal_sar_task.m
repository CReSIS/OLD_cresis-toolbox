function success = rx_chan_equal_sar_task(param)
% success = rx_chan_equal_sar_task(param)
%
% Task/job running on cluster that is called from rx_chan_equal_sar. This
% function is generally not called directly.
%
% Inputs:
% data = Nt by Nx by Nc single/double matrix
%   Nt = fast time
%   Nx = slow time
%   Nc = cross-track channels
% param = structure controlling receiver channel equalization
%  .rx_paths = single wf-adc list (Nc by 2 matrix)
%    List of wf-adc pairs.
%    [wf adc; wf adc; ... ; wf adc];
%  .td = time delay correction for each wf_adc pair (Nc by 1 vector)
%    e.g. -12 ns implies channel is delayed by 12 ns so 12 ns of delay
%    will be removed (targets will move to a closer range)
%    default is all zeros, units of seconds
%  .phase = phase correction for each wf_adc pair (Nc by 1 vector)
%    e.g. 35 deg implies channel leads by 35 deg so -35 deg phase will be
%    applied (target phases will have -35 deg phase from original)
%    default is all zeros, units of degrees (angle(voltage)*180/pi)
%  .amp = amplitude correction for each wf_adc pair (Nc by 1 vector)
%    e.g. 2 dB implies channel is 2 dB larger so -2 dB amplitude adjustment
%    will be applied (targets will be 2 dB smaller)
%    default is all zeros, units of log power (20*log(voltage))
%  .ref_bins = bin range around the target to use in cross correlation
%    [-12 12] searches 12 bins before and after
%    default is [-12 12]
%  .search_bins = bin range to search for corresponding target in other
%    channels, [-7 7] searches 7 bins forward and backward allowing the
%    channels to have time delay offsets up to 7 range bins
%    default is [-7 7]
%  .Mt = over-sampling factor to apply in time domain
%    default is 100
%  .plot_en = enable plotting of outputs
%    default is false
%  .cross_correlation_flag = enable use of cross correlation instead of
%    just peak finding (required for obtaining td values)
%  .freq = double Nt by 1 vector, fast-time axis of data, Hz
%  .time = double Nt by 1 vector, fast-time axis of data, seconds
%
%
% Outputs:
%   td_out = recommended hdr.td
%   phase_out = recommended hdr.phase
%   amp_out = recommended hdr.amp
%
% Authors: Peng Seng Tan, John Paden

physical_constants;

fprintf('=====================================================================\n');
fprintf('Rx Equalization Task %s_%03d chunks %d to %d (%s)\n\n', ...
  param.day_seg, param.equal.frm, param.equal.chunks(1), ...
  param.equal.chunks(end), datestr(now,'HH:MM:SS'));

param.equal.Nsv = 64;
param.equal.theta_rng = [-10 10]/180*pi;

% ======================================================================
% Input arguments and general setup
% ======================================================================

% =====================================================================
% Load the data
% =====================================================================

if strcmpi(param.csarp.sar_type,'f-k')
  sar_type = 'fk';
else
  sar_type = 'tdc';
end
in_path = fullfile(ct_filename_out(param, ...
  param.equal.in_path, 'CSARP_out'), ...
  sprintf('%s_data_%03d_01_01',sar_type,param.equal.frm));

chunk = param.equal.chunks(1);
data = [];
chan_equal = [];
fcs = [];
for img = 1:length(param.equal.imgs)
  
  wf_adc_list = param.equal.imgs{img};
  for wf_adc_idx = 1:size(wf_adc_list,1)
    wf = wf_adc_list(wf_adc_idx,1);
    adc = wf_adc_list(wf_adc_idx,2);
    
    fn = fullfile(in_path,sprintf('wf_%02d_adc_%02d_chk_%03d.mat', ...
      wf, adc, chunk));
    if param.debug_level >= 1
      fprintf('  Loading %s (%s)\n', fn, datestr(now,'HH:MM:SS'));
    end
    sar_data = load(fn);
    chan_equal{img}(wf_adc_idx) = 10.^(param.radar.wfs(wf).chan_equal_dB(param.radar.wfs(wf).rx_paths(adc))/20) ...
      .* exp(j*param.radar.wfs(wf).chan_equal_deg(param.radar.wfs(wf).rx_paths(adc))/180*pi);
    data{img}(:,1:size(sar_data.fk_data,2),wf_adc_idx) = sar_data.fk_data;
    fcs{img}{wf_adc_idx} = sar_data.fcs;
    wfs = sar_data.wfs;
  end
end

%% Contact ops to get surface and bottom

if ~isfield(param.combine,'layer_fn') || isempty(param.combine.layer_fn)
  param.combine.layer_fn = 'layerData';
end

%% Get the generic layer data path
layer_path = fullfile(ct_filename_out(param,param.combine.layer_fn,'',0));

%% Load the current frame
layer_fn = fullfile(layer_path,sprintf('Data_%s_%03d.mat',param.day_seg,param.equal.frm));
layer = load(layer_fn);
layer_gps_time = layer.GPS_time;
surface = layer.layerData{1}.value{2}.data;
bottom = layer.layerData{2}.value{2}.data;
%% Get the previous frame if necessary
if fcs{1}{1}.gps_time(1) < layer_gps_time(1)-1
  layer_fn = fullfile(layer_path,sprintf('Data_%s_%03d.mat',param.day_seg,param.equal.frm-1));
  if exist(layer_fn,'file')
    layer = load(layer_fn);
    layer_gps_time = [layer.GPS_time layer_gps_time];
    surface = [layer.layerData{1}.value{2}.data surface];
    bottom = [layer.layerData{1}.value{2}.data bottom];
  end
end
%% Get the next frame if necessary
if fcs{1}{1}.gps_time(end) > layer_gps_time(end)+1
  layer_fn = fullfile(layer_path,sprintf('Data_%s_%03d.mat',param.day_seg,param.equal.frm+1));
  if exist(layer_fn,'file')
    layer = load(layer_fn);
    layer_gps_time = [layer_gps_time layer.GPS_time];
    surface = [surface layer.layerData{1}.value{2}.data];
    bottom = [bottom layer.layerData{1}.value{2}.data];
  end
end
%% Since layer files may have overlapping data, sort it
[layer_gps_time new_surface_idxs] = sort(layer_gps_time);
surface = surface(new_surface_idxs);
bottom = bottom(new_surface_idxs);

%% Do the interpolation
good_idxs = find(isfinite(surface));
surface = interp1(layer_gps_time(good_idxs),surface(good_idxs),fcs{1}{1}.gps_time,'linear','extrap');
good_idxs = find(isfinite(bottom));
if length(good_idxs) < 2
  bottom = NaN*zeros(size(fcs{1}{1}.gps_time));
else
  bottom = interp1(layer_gps_time(good_idxs),bottom(good_idxs),fcs{1}{1}.gps_time,'linear','extrap');
end
bottom = surface;

bottom_bins = round(interp1(wfs(wf).time, ...
  1:length(wfs(wf).time), bottom));

% Interpolate bed data onto range bins
param.equal.layer(1,:) = round(interp1(wfs(wf).time, ...
  1:length(wfs(wf).time), bottom)) - 6;
param.equal.layer(1,isnan(param.equal.layer(1,:))) = 1;

param.equal.layer(2,:) = round(interp1(wfs(wf).time, ...
  1:length(wfs(wf).time), bottom)) + 6;
param.equal.layer(2,~isnan(param.equal.layer(2,:))) ...
  = param.equal.layer(1,~isnan(param.equal.layer(2,:))) ...
  + param.equal.layer(2,~isnan(param.equal.layer(2,1))) ...
  - param.equal.layer(1,~isnan(param.equal.layer(1,1)));
param.equal.layer(2,isnan(param.equal.layer(2,:))) = length(wfs(wf).time);

for img = 1:length(param.equal.imgs)
  wf_adc_list = param.equal.imgs{img};
  wf_adc_idx = 1;
  wf = wf_adc_list(wf_adc_idx,1);
  adc = wf_adc_list(wf_adc_idx,2);
  
  Time = sar_data.wfs(wf).time;
  %fcs{img}{wf_adc_idx}.bottom = fcs{img}{wf_adc_idx}.surface + 34e-6;
  
  % figure(1); clf;
  % imagesc([],Time,angle(data{img}(:,:,1) .* conj(data{img}(:,:,3))))
  % hold on;
  % plot(fcs{img}{wf_adc_idx}.surface,'k')
  % plot(fcs{img}{wf_adc_idx}.bottom,'k')
  % hold off;
  
  % Number of fast-time samples in the data
  Nt = size(data{img},1);
  
  % Number of slow-time/along-track samples in the data
  Nx = size(data{img},2);
  
  % Number of cross-track channels in the data
  Nc = size(data{img},3);
  
  %% Apply channel equalization coefficients.
  for chan = 1:Nc
    data{img}(:,:,chan) = data{img}(:,:,chan) / chan_equal{img}(chan);
  end
  
  %% Refine layer
  surf_data = fir_dec(abs(mean(data{img},3)).^2,ones(1,9),1);
  new_bottom_bins = bottom_bins;
  for rline=1:size(surf_data,2)
    [~,offset] = max(surf_data(bottom_bins(rline)+(-5:5),rline));
    new_bottom_bins(rline) = bottom_bins(rline) - 6 + offset;
  end
  
%   imagesc(lp(surf_data));
%   hold on
%   plot(bottom_bins,'k');
%   plot(new_bottom_bins,'g');
%   hold off;
  
  bottom_bins = new_bottom_bins;
  
  %% Set up which range bins will be output bins
%   if img == 0
%     bins = 300 ...
%       : param.equal.dbin ...
%       : 550;
%   elseif img == 1
%         bins = 800 ...
%           : param.equal.dbin ...
%           : 900;
% %     bins = 700 ...
% %       : param.equal.dbin ...
% %       : 800;
%     %     bins = 1330 ...
%     %       : param.equal.dbin ...
%     %       : 1400;
%   else
%     bins = numel(param.equal.bin_rng)/2+0.5 ...
%       : param.equal.dbin ...
%       : Nt-(numel(param.equal.bin_rng)/2-0.5);
%   end
  bins = min(param.equal.layer(1,:)) : param.equal.dbin : max(param.equal.layer(2,:));

  %% Set up which range lines will be output lines
  
  % Start/stop output range lines passed in (typical operation from
  % combine_wf_chan_task)
  rlines = param.equal.rlines(1): param.equal.dline ...
    : min(param.equal.rlines(2),size(data{img},2)-max(param.equal.rline_rng));
  
  %% Preallocate output matrices
  % equal_Rxx = complex(zeros(Nc^2,length(bins),length(rlines)));
  equal_pow = zeros(Nc,length(bins),length(rlines));
  equal_angle = complex(zeros(Nc,length(bins),length(rlines)));
  Smusic = zeros(param.equal.Nsv,length(bins),length(rlines));
  Smvdr = zeros(param.equal.Nsv,length(bins),length(rlines));
  
  %% Setup Steering Vector variables
  theta = param.equal.sv_fh(param.equal.Nsv,wfs(1).fc);
  % Beam-forming approach: Determine which ky indices will be used in
  % selecting the value for each range bin
  param.equal.freq_rng = find(theta >= param.equal.theta_rng(1) & theta <= param.equal.theta_rng(2));
  if isempty(param.equal.freq_rng)
    [tmp param.equal.freq_rng] = min(abs(theta-mean(param.equal.theta_rng)));
  end
  
  %% Estimate equalization at each pixel by looping through each range line and each
  % range bin
  tic;
  rline = 0;
  for rline_idx = 1:1:length(rlines)
    if param.debug_level >= 3
      % Insert debug code as needed
      if rline_idx == inf
        figure(1); clf;
        imagesc(Smusic);
        figure(2); clf;
        imagesc(20*log10(abs(squeeze(data{img}(:,line,:))).^2));
        keyboard
      end
    end
    
    rline = rlines(rline_idx);
    if param.debug_level >= 2
      fprintf('    Record %.0f (%.0f of %.0f) (%.1f sec)\n', rline, rline_idx, ...
        length(rlines), toc);
    end
    
    %% Check for edge conditions (not enough data to create dataSample
    % properly)
    if rline+param.equal.rline_rng(1) < 1
      rline_rng = 1-rline : param.equal.rline_rng(end);
    elseif  rline+param.equal.rline_rng(end) > Nx
      rline_rng = param.equal.rline_rng(1) : Nx-rline;
    else
      rline_rng = param.equal.rline_rng;
    end
    
    %% Steering Vector setup
    % Make column vectors of y and z-positions
    if 1
      for wf_adc_idx = 1:length(fcs{img})
        y_pos(wf_adc_idx,1) = fcs{img}{wf_adc_idx}.pos(2,rline);
        z_pos(wf_adc_idx,1) = fcs{img}{wf_adc_idx}.pos(3,rline);
      end
    else
      % FOR DEBUGGING... assumes level flight
      for wf_adc_idx = 1:length(fcs{img})
        wf = wf_adc_list(wf_adc_idx,1);
        adc = wf_adc_list(wf_adc_idx,2);
        param.gps_source = fcs{1}{1}.gps_source;
        phase_center = lever_arm(param, param.radar.wfs(wf).tx_weights, param.radar.wfs(wf).rx_paths(adc));
        y_pos(wf_adc_idx,1) = -phase_center(2);
        z_pos(wf_adc_idx,1) = -phase_center(3);
      end
    end
    % Determine Steering Vector for beam-forming approaches
    % Make column vectors of y and z-positions
    [theta param.equal.sv{img}] = param.equal.sv_fh(param.equal.Nsv,wfs(wf).fc,y_pos,z_pos);
    
    %% Iterate through each range bin
    for bin_idx = 1:1:length(bins)
      bin = bins(bin_idx);
      
      %% Get data sample (snapshots/multilooks) LxM L=#looks M=#sensors
      dataSample = data{img}(bin+param.equal.bin_rng,rline+rline_rng,:);
      dataSample = reshape(dataSample,[length(param.equal.bin_rng)*length(rline_rng) Nc]).';
      
      %% Estimate direction of arrival
      Rxx = 1/size(dataSample,1) * (dataSample * dataSample');
      [V,D] = eig(Rxx);
      eigenVals = diag(D);
      [eigenVals noiseIdxs] = sort(eigenVals);
      param.equal.Nsig = 1;
      noiseIdxs = noiseIdxs(1:end-param.equal.Nsig);
      Smusic(:,bin_idx,rline_idx) = fftshift(1./mean(abs(param.equal.sv{img}(:,:)'*V(:,noiseIdxs)).^2,2),1);
      diagonal = sqrt(mean(mean(abs(Rxx).^2))) * diag(ones(Nc,1),0);
      param.equal.diag_load = 0;
      warning off;
      Rxx_inv = inv(Rxx + param.equal.diag_load*diagonal);
      warning on;
      for freq_idx = 1:size(param.equal.sv{1},2)
        Smvdr(freq_idx,bin_idx,rline_idx) = 1 ./ real(param.equal.sv{img}(:,freq_idx)' * Rxx_inv * param.equal.sv{img}(:,freq_idx));
      end
      Smvdr(:,bin_idx,rline_idx) = fftshift(Smvdr(:,bin_idx,rline_idx));
      
      %% Estimate power and relative phase/amplitude of data samples AFTER
      % applying steering vector correction for nadir direction
      dataSample = dataSample .* repmat(conj(param.equal.sv{img}(:,1)),[1 size(dataSample,2)]);
      %equal_Rxx(:,bin_idx,rline_idx) = Rxx(:);
      equal_pow(:,bin_idx,rline_idx) = dot(dataSample.',dataSample.');
      equal_angle(:,bin_idx,rline_idx) = dot(repmat(dataSample(param.equal.ref_wf_adc_idx,:),[Nc 1]).', dataSample.');
      
            if 0 && img == 1 && rline_idx >= 10 && bin_idx == length(bins) %% DEBUGGING
              figure(1); clf;
              plot(fftshift(theta)*180/pi, Smusic(:,bin_idx,rline_idx));
              hold on;
              plot(fftshift(theta)*180/pi, lp(Smvdr(:,bin_idx,rline_idx)) - max(lp(Smvdr(:,bin_idx,rline_idx))) ,'r');
              hold off;
              equal_angle(:,bin_idx,rline_idx) = dot(repmat(dataSample(:,param.equal.ref_wf_adc_idx), [1 Nc]), dataSample);
              grid on;
              figure(2); clf;
%               plot(angle(equal_angle(:,bin_idx,rline_idx)).'*180/pi,'r')
%               hold on;
%               plot(angle(param.equal.sv{img}(:,1)).'*180/pi,'k')
%               hold off;
              bottom_bin_idx = interp1(bins, 1:length(bins), bottom_bins(rline),'nearest');
              imagesc(angle(equal_angle(:,:,rline_idx)))
              hold on;
              plot([bottom_bin_idx bottom_bin_idx], [1 Nc],'k');
              hold off;
              keyboard
              % imagesc(angle(data{img}(:,:,7) .* conj(data{img}(:,:,param.equal.ref_wf_adc_idx)) ))
            end
    end
  end
  
  % figure(1); clf;
  % imagesc(lp(squeeze(equal_angle(1,:,:))))
  % colorbar;
  % figure(2); clf;
  % imagesc(angle(squeeze(equal_angle(1,:,:))))
  % colorbar;
  %
  % set(2,'Position',get(1,'Position'));
  
  equal_angle_layer = zeros(size(equal_angle,1),size(equal_angle,3));
  % For each range line find the closest bin to the layer
  bottom_bin_idx = interp1(bins, 1:length(bins), bottom_bins,'nearest');
  for rline = 1:size(equal_angle,3)
    equal_angle_layer(:,rline) = equal_angle(:,bottom_bin_idx(rline),rline);
  end
  
  chan_equal_est = equal_angle_layer;
  fprintf('%.1f\t', angle(mean(chan_equal_est,2)) * 180/pi);
  fprintf('\n');
  
  theta_fftshift = fftshift(theta);
  
  roll_est = zeros(1,size(Smusic,3));
  good_theta_idxs = find(abs(theta_fftshift) < 20/180*pi);
  for rline = 1:size(Smusic,3)
    [music_peak_val music_peak_idx] = max(medfilt2(Smusic(good_theta_idxs,:,rline),[3 5]));
    music_peak_idx = music_peak_idx + good_theta_idxs(1)-1;
    roll_est(rline) = median(theta_fftshift(music_peak_idx));
  end
  [B,A] = butter(2,0.2);
  roll_est_filt = filtfilt(B,A,roll_est);
  figure(1); clf;
  %   plot(fcs{img}{1}.gps_time(rlines), (roll_est_filt - fcs{img}{1}.roll(rlines))*180/pi)
  plot(fcs{img}{1}.gps_time(rlines), roll_est_filt*180/pi)
  hold on;
  plot(fcs{img}{1}.gps_time, fcs{img}{1}.roll*180/pi,'r');
  hold off;
  legend('Estimated','Actual');
  xlabel('GPS time (sec)');
  ylabel('Roll (deg)');
  
  drawnow;
  
  out_fn = fullfile(ct_filename_out(param,param.equal.out_path,'CSARP_out'), ...
    sprintf('equal_%03d', param.equal.frm), sprintf('chk_%03d_img_%02d.mat',param.equal.chunks,img));
  out_fn_dir = fileparts(out_fn);
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  param_equal = param;
  param_equal.wfs = sar_data.wfs;
  param_records = sar_data.param_records;
  param_csarp = sar_data.param_csarp;
  
  [Latitude,Longitude,Elevation] ...
    = ecef2geodetic(sar_data.fcs.origin(1,rlines), ...
    sar_data.fcs.origin(2,rlines), ...
    sar_data.fcs.origin(3,rlines),WGS84.ellipsoid);
  Latitude = Latitude*180/pi;
  Longitude = Longitude*180/pi;
  GPS_time = sar_data.fcs.gps_time(rlines);
  Surface = sar_data.fcs.surface(rlines);
  Bottom = sar_data.fcs.bottom(rlines);
  Roll = sar_data.fcs.roll(rlines);
  Pitch = sar_data.fcs.pitch(rlines);
  Heading = sar_data.fcs.heading(rlines);
  
  fprintf('  Writing %s\n', out_fn);
  save(out_fn, 'Smusic', 'equal_pow', 'equal_angle', 'roll_est', ...
    'chan_equal_est', 'GPS_time', 'Latitude', 'Longitude', 'Elevation', ...
    'Roll', 'Pitch', 'Heading', 'Surface', 'Bottom', 'param_equal', ...
    'param_records', 'param_csarp');
end

success = true;

return;

